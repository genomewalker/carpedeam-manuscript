library(tidyverse)
library(janitor)
library(showtext)
library(tidygraph)
library(igraph)
library(UpSetR)

showtext_auto()

# Get contigs stats
contigs_files <- list.files("./data/genes/ancientGut", full.names = TRUE, pattern = ".contig_stats.tsv")

read_contig_files <- function(filename) {
    pattern <- "^(.*)_(.*)\\.contig_stats\\.tsv$"
    matches <- str_match(basename(filename), pattern)

    # Extract the matched groups
    sample <- matches[2]
    assembler <- matches[3]
    df <- read_tsv(filename, col_names = FALSE) |>
        setNames(c("contig_id", "length")) |>
        mutate(sample = sample, assembler = assembler)
    if (nrow(df) == 0) {
        return(NULL)
    }
    return(df)
}

contigs <- map_dfr(contigs_files, read_contig_files)

min_contig_length <- 1000

contig_stats <- contigs |>
    filter(length >= min_contig_length) |>
    group_by(sample, assembler) |>
    summarise(
        n_contigs = n(),
        total_length = sum(length),
        mean_length = mean(length),
        median_length = median(length),
        min_length = min(length),
        max_length = max(length),
        .groups = "drop"
    )

# HOw similar are the different genomes used for simulating the ancient gut?
skani_allvsall_genome <- read_tsv("./data/genomes/ancientGut/gut_simulation_allvsall.tsv") |>
    clean_names() |>
    rowwise() |>
    mutate(avg_af = ((align_fraction_ref + align_fraction_query) / 2 / 100)) |>
    arrange(desc(ani), desc(avg_af))

skani_allvsall_genome |>
    filter(ani >= 95) |>
    ggplot(aes(x = avg_af, y = "af")) +
    geom_boxplot() +
    geom_jitter(aes(fill = ani / 100), width = 0.1, size = 3, shape = 21) +
    theme_bw() +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1, "cm")
    ) +
    xlab("Average alignment fraction") +
    ylab("") +
    scale_x_continuous(labels = scales::percent) +
    scale_fill_viridis_c(name = "ANI", labels = scales::percent_format(accuracy = 1)) +
    coord_fixed(ratio = 0.03)


# Which is the level of redundancy in the different assemblies based on the ANI?
skani_contigs_files <- list.files("./data/genes/ancientGut", full.names = TRUE, pattern = "\\.ani.tsv")

read_skani_contig_file <- function(filename) {
    pattern <- "^(.*)_(.*)\\.ani\\.tsv$"
    matches <- str_match(basename(filename), pattern)

    # Extract the matched groups
    sample <- matches[2]
    assembler <- matches[3]
    df <- read_tsv(filename) |>
        clean_names() |>
        mutate(sample = sample, assembler = assembler)
    if (nrow(df) == 0) {
        return(NULL)
    }
    return(df)
}

skani_contig <- map_dfr(skani_contigs_files, read_skani_contig_file)

skani_contig_edges <- skani_contig |>
    inner_join(contigs |> rename(ref_name = contig_id, ref_length = length)) |>
    inner_join(contigs |> rename(query_name = contig_id, query_length = length)) |>
    filter(query_length >= min_contig_length, ref_length >= min_contig_length) |>
    filter(ani >= 99) |>
    filter(if_else(ref_length <= query_length, align_fraction_ref, align_fraction_query) >= 90) |>
    select(sample, assembler, ref_name, query_name, ani, align_fraction_ref, align_fraction_query, ref_length, query_length)

s_a <- skani_contig_edges |>
    select(sample, assembler) |>
    distinct()


# Find non-overlapping cliques
find_non_overlapping_cliques <- function(cliques) {
    non_overlapping_cliques <- list()
    used_nodes <- c()

    for (clique in cliques) {
        if (length(intersect(clique, used_nodes)) == 0) {
            non_overlapping_cliques <- append(non_overlapping_cliques, list(clique))
            used_nodes <- union(used_nodes, clique)
        }
    }

    return(non_overlapping_cliques)
}

# Create a graph for each sample and assembler
g <- s_a |>
    mutate(graph = map2(sample, assembler, ~ {
        # Create the graph
        graph <- skani_contig_edges |>
            filter(sample == .x, assembler == .y) |>
            select(ref_name, query_name, weight = ani, align_fraction_ref, align_fraction_query) |>
            as_tbl_graph(directed = FALSE) |>
            mutate(component = group_components()) |>
            inner_join(contigs |> filter(sample == .x, assembler == .y) |> select(name = contig_id, length))

        # Detect cliques and find non-overlapping cliques
        cliques <- max_cliques(graph)
        non_overlapping_cliques <- find_non_overlapping_cliques(cliques)

        # Assign clique membership to nodes
        membership <- rep(NA, vcount(graph))
        for (i in seq_along(non_overlapping_cliques)) {
            for (node in non_overlapping_cliques[[i]]) {
                membership[node] <- i
            }
        }

        # Add clique membership to the graph
        graph <- graph |>
            activate(nodes) |>
            mutate(clique_membership = membership)

        return(graph)
    }))


# Function to get components and their sizes
get_component_info <- function(graph) {
    components <- graph |>
        as_tibble() |>
        group_by(component) |>
        summarise(size = n(), nodes = list(name)) |>
        arrange(desc(size))

    num_components <- n_distinct(components$component)

    list(components = components, num_components = num_components, num_nodes = igraph::vcount(graph))
}

# Apply the function to each graph in the list and get stats
graph_stats <- g |>
    mutate(
        component_info = map(graph, get_component_info),
        component_size = map(component_info, "components"),
        num_components = map_int(component_info, "num_components"),
        num_nodes = map_int(component_info, "num_nodes"),
    ) |>
    inner_join(contig_stats |> select(sample, assembler, n_contigs), by = c("sample", "assembler"))

# Where do this contigs map in the reference genomes?
skani_search_files <- list.files("./data/genes/ancientGut", full.names = TRUE, pattern = "\\.search-skani.tsv")
read_skani_search_file <- function(filename) {
    pattern <- "^(.*)_(.*)\\.search-skani\\.tsv$"
    matches <- str_match(basename(filename), pattern)

    # Extract the matched groups
    sample <- matches[2]
    assembler <- matches[3]
    df <- read_tsv(filename) |>
        clean_names() |>
        mutate(sample = sample, assembler = assembler)
    if (nrow(df) == 0) {
        return(NULL)
    }
    return(df)
}

# Get the data of mapping the contigs to the original assemblies
skani_search <- map_dfr(skani_search_files, read_skani_search_file) |>
    filter(align_fraction_query >= 95, ani >= 99) |>
    inner_join(contigs |> rename(query_name = contig_id, qlen = length), by = c("sample", "assembler", "query_name")) |>
    filter(qlen > min_contig_length)


# Combine search results with grpah nodes
join_nodes_with_search <- function(graph, sample, assembler) {
    hits <- skani_search |>
        filter(sample == sample, assembler == assembler) |>
        group_by(query_name) |>
        filter(
            ani == max(ani) &
                align_fraction_query == max(align_fraction_query)
        ) |>
        slice(1) |>
        ungroup() |>
        select(query_name, ref_name) |>
        distinct()


    graph |>
        activate(nodes) |>
        as_tibble() |>
        rename(query_name = name) |>
        left_join(hits) |>
        mutate(is_mapped = ifelse(!is.na(ref_name), "mapped", "unmapped")) |>
        # filter(component == 1, clique_membership == 28) |>
        arrange(clique_membership, desc(length)) |>
        group_by(ref_name, component, clique_membership) |>
        summarise(
            contig_num = n(),
            max_contig_bp = max(length),
            sum_contig_bp = sum(length),
            unmapped_bp = sum(ifelse(is.na(ref_name), length, 0)),
            .groups = "drop"
        ) |>
        mutate(is_mapped = ifelse(!is.na(ref_name), "mapped", "unmapped")) |>
        group_by(clique_membership) |>
        add_count(name = "ref_num") |>
        ungroup() |>
        mutate(type = ifelse(ref_num > 1, "multi", "single"))
}

# Apply the function to each row in graph_stats
graph_stats <- graph_stats |>
    mutate(joined_nodes = pmap(list(graph, sample, assembler), join_nodes_with_search))


# Calculate for each assembly|sample the number of contigs that map to more than one reference genome
# calculate the number of contigs that are duplicated and need to be removed
get_ncontigs_multiple_refs <- function(joined_nodes) {
    df <- joined_nodes |>
        mutate(type = ifelse(is_mapped == "unmapped", "unmapped", type)) |>
        group_by(type) |>
        summarise(
            n = n(),
            contig_num = sum(contig_num),
            contig_unmmaped_bp = sum(unmapped_bp),
            contig_mapped_bp = sum(sum_contig_bp) - contig_unmmaped_bp,
            contig_bp_total = sum(sum_contig_bp),
            contig_bp_keep = sum(max_contig_bp),
            contig_bp_rm = contig_bp_total - contig_bp_keep,
            .groups = "drop"
        ) |>
        mutate(contig2rm = ifelse(type == "multi", 0, contig_num - n))
}

graph_stats <- graph_stats |>
    mutate(
        node_summary = map(joined_nodes, get_ncontigs_multiple_refs),
        n_contigs_multi = map_int(node_summary, ~ {
            # Filter rows where type is "multi"
            multi_rows <- .x |> filter(type == "multi")
            # Extract the first contig_num of the filtered rows (if any)
            if (nrow(multi_rows) > 0) multi_rows$contig_num[1] else 0
        }),
        n_contigs_single = map_int(node_summary, ~ {
            # Filter rows where type is "single"
            single_rows <- .x |> filter(type == "single")
            # Extract the first contig_num of the filtered rows (if any)
            if (nrow(single_rows) > 0) single_rows$contig_num[1] else 0
        }),
        n_contigs_unmapped = map_int(node_summary, ~ {
            # Filter rows where type is "single"
            single_rows <- .x |> filter(type == "unmapped")
            # Extract the first contig_num of the filtered rows (if any)
            if (nrow(single_rows) > 0) single_rows$contig_num[1] else 0
        }),
        n_contigs_single_rm = map_int(node_summary, ~ {
            # Filter rows where type is "single"
            single_rows <- .x |> filter(type == "single")
            # Extract the first contig_num of the filtered rows (if any)
            if (nrow(single_rows) > 0) single_rows$contig2rm[1] else 0
        }),
        n_contigs_multi_rm = map_int(node_summary, ~ {
            # Filter rows where type is "multi"
            single_rows <- .x |> filter(type == "multi")
            # Extract the first contig_num of the filtered rows (if any)
            if (nrow(single_rows) > 0) single_rows$contig2rm[1] else 0
        }),
        n_contigs_unmapped_rm = map_int(node_summary, ~ {
            # Filter rows where type is "single"
            single_rows <- .x |> filter(type == "unmapped")
            # Extract the first contig_num of the filtered rows (if any)
            if (nrow(single_rows) > 0) single_rows$contig2rm[1] else 0
        }),
        contig_bp_single_keep = map_int(node_summary, ~ {
            # Filter rows where type is "single"
            single_rows <- .x |> filter(type == "single")
            # Extract the first contig_num of the filtered rows (if any)
            if (nrow(single_rows) > 0) single_rows$contig_bp_keep[1] else 0
        }),
        contig_bp_single_rm = map_int(node_summary, ~ {
            # Filter rows where type is "single"
            single_rows <- .x |> filter(type == "single")
            # Extract the first contig_num of the filtered rows (if any)
            if (nrow(single_rows) > 0) single_rows$contig_bp_rm[1] else 0
        }),
        contig_bp_multi_keep = map_int(node_summary, ~ {
            # Filter rows where type is "multi"
            single_rows <- .x |> filter(type == "multi")
            # Extract the first contig_num of the filtered rows (if any)
            if (nrow(single_rows) > 0) single_rows$contig_bp_keep[1] else 0
        }),
        contig_bp_multi_rm = map_int(node_summary, ~ {
            # Filter rows where type is "multi"
            single_rows <- .x |> filter(type == "multi")
            # Extract the first contig_num of the filtered rows (if any)
            if (nrow(single_rows) > 0) single_rows$contig_bp_rm[1] else 0
        }),
        contig_bp_unmapped_keep = map_int(node_summary, ~ {
            # Filter rows where type is "single"
            single_rows <- .x |> filter(type == "unmapped")
            # Extract the first contig_num of the filtered rows (if any)
            if (nrow(single_rows) > 0) single_rows$contig_bp_keep[1] else 0
        }),
        contig_bp_unmapped_rm = map_int(node_summary, ~ {
            # Filter rows where type is "single"
            single_rows <- .x |> filter(type == "unmapped")
            # Extract the first contig_num of the filtered rows (if any)
            if (nrow(single_rows) > 0) single_rows$contig_bp_rm[1] else 0
        })
    )

# Get all stats
mapping_stats <- skani_search |>
    select(sample, assembler, query_name, qlen) |>
    distinct() |>
    group_by(sample, assembler) |>
    summarise(mapped_contigs = n(), mapped_bp = sum(qlen), .groups = "drop") |>
    inner_join(contig_stats |> select(sample, assembler, n_contigs, total_length), by = c("sample", "assembler")) |>
    left_join(graph_stats |>
        select(
            sample,
            assembler,
            n_contigs_multi,
            n_contigs_single,
            n_contigs_unmapped,
            n_contigs_single_rm,
            n_contigs_multi_rm,
            n_contigs_unmapped_rm,
            contig_bp_single_keep,
            contig_bp_single_rm,
            contig_bp_multi_keep,
            contig_bp_multi_rm,
            contig_bp_unmapped_keep,
            contig_bp_unmapped_rm
        ), by = c("sample", "assembler")) |>
    replace_na(list(
        n_contigs = 0,
        n_contigs_multi = 0,
        n_contigs_single = 0,
        n_contigs_unmapped = 0,
        n_contigs_single_rm = 0,
        n_contigs_multi_rm = 0,
        n_contigs_unmapped_rm = 0,
        contig_bp_single_keep = 0,
        contig_bp_single_rm = 0,
        contig_bp_multi_keep = 0,
        contig_bp_multi_rm = 0,
        contig_bp_unmapped_keep = 0,
        contig_bp_unmapped_rm = 0
    ))


# Do figure 2
colors <- c(
    "mapped_nondup_bp" = "#2066a8",
    "mapped_dup_var_bp" = "#3594cC",
    "mapped_dup_red_rep_bp" = "#8cc5e3",
    "mapped_dup_red_bp" = "#F09868FF",
    "unmapped_dup_bp" = "#c46666",
    "unmapped_nondup_bp" = "#d8a6a6"
)

mapping_stats |>
    mutate(
        # Total mapped bp
        mapped_bp = mapped_bp,
        # Total unmapped bp
        unmapped_bp = total_length - mapped_bp,
        # Total unmapped bp that are not duplicates
        unmapped_nondup_bp = unmapped_bp - (contig_bp_unmapped_keep + contig_bp_unmapped_rm),
        # Total unmapped bp that are duplicates
        unmapped_dup_bp = unmapped_bp - unmapped_nondup_bp,
        # Total mapped bp that are duplicates and can be due to variation
        mapped_dup_var_bp = contig_bp_multi_keep,
        # Total mapped bp that are duplicates and can be due to redundancy
        mapped_dup_red_bp = contig_bp_single_rm + contig_bp_multi_rm,
        # Total mapped bp that from the duplicates that are representative
        mapped_dup_red_rep_bp = contig_bp_single_keep,
        # Total mapped bp that are duplicates
        mapped_dup_bp = mapped_dup_var_bp + mapped_dup_red_bp + mapped_dup_red_rep_bp,
        # Total mapped bp that are not duplicates
        mapped_nondup_bp = mapped_bp - mapped_dup_bp,
    ) |>
    select(
        # total_length,
        # mapped_bp,
        # unmapped_bp,
        sample,
        assembler,
        mapped_nondup_bp,
        mapped_dup_var_bp,
        mapped_dup_red_bp,
        mapped_dup_red_rep_bp,
        unmapped_nondup_bp,
        unmapped_dup_bp
    ) |>
    # mutate(
    #     mapped_sub = mapped_dup_var_bp + mapped_dup_red_bp + mapped_dup_red_rep_bp + mapped_nondup_bp,
    #     unmapped_sub = unmapped_dup_bp + unmapped_nondup_bp,
    #     total_sub = mapped_sub + unmapped_sub
    # ) |>
    # select(sample, assembler, total_length, total_sub, mapped_bp, mapped_sub, unmapped_bp, unmapped_sub)
    pivot_longer(cols = c(
        mapped_dup_var_bp,
        mapped_dup_red_rep_bp,
        mapped_dup_red_bp,
        mapped_nondup_bp,
        unmapped_dup_bp,
        unmapped_nondup_bp
    ), names_to = "name", values_to = "bp") |>
    mutate(type = ifelse(grepl("mapped", name), "mapped", "unmapped")) |>
    mutate(
        cov = case_when(
            grepl("c3", sample) ~ "3X",
            grepl("c5", sample) ~ "5X",
            grepl("c10", sample) ~ "10X",
        ),
        dmg = case_when(
            grepl("mid", sample) ~ "mid",
            grepl("high", sample) ~ "high",
        ),
        assm = case_when(
            grepl("carpedeam-unsafe", assembler) ~ "Carpedeam\n(unsafe)",
            grepl("carpedeam-safe", assembler) ~ "Carpedeam\n(safe)",
            grepl("megahit", assembler) ~ "Megahit",
        )
    ) |>
    mutate(
        name = fct_relevel(name, c("unmapped_dup_bp", "unmapped_nondup_bp", "mapped_dup_red_bp", "mapped_dup_red_rep_bp", "mapped_dup_var_bp", "mapped_nondup_bp")),
        cov = fct_relevel(cov, c("3X", "5X", "10X")),
        dmg = fct_relevel(dmg, c("high", "mid"))
    ) |>
    ggplot(aes(x = assm, y = bp, fill = name)) +
    geom_bar(stat = "identity", position = "stack", width = 1, color = "black", linewidth = 0.3) +
    scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6)) +
    scale_fill_manual(values = colors, name = NULL) +
    facet_grid2(c("dmg", "cov"),
        labeller = "label_both",
        scales = "free_y", independent = "y"
    ) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        legend.position = "bottom",
    ) +
    ylab("Base pairs") +
    xlab("")


### PROTEIN analyses
# How many proteins are in the different assemblies?
files <- list.files("./data/genes/ancientGut", full.names = TRUE, pattern = "protein_stats")
read_protein_stats <- function(X) {
    df <- read_tsv(X) |>
        clean_names()
    return(df)
}

# Get all the protein stats
protein_stats <- map_dfr(files, read_protein_stats) |>
    mutate(identity = ifelse(is.na(identity), "raw", identity)) |>
    mutate(identity = ifelse(identity == "1", "1.0", identity)) |>
    select(-file, -format) |>
    separate(sample, into = c("sample", "assembler"), sep = "_(?!.*_)", extra = "merge")

# Let's plot the amount of reduction from raw to dereplicated. It should be
data <- protein_stats %>%
    select(assembler, sample, identity, num_seqs) %>%
    pivot_wider(names_from = identity, values_from = num_seqs, values_fill = 0) |>
    mutate(
        c95 = `0.95` / raw,
        c1 = (`1.0` / raw) - c95,
        rawp = 1 - (c95 + c1),
        c95 = c95 * raw,
        c1 = c1 * raw,
        rawp = raw * rawp
    ) |>
    select(assembler, sample, c95, c1, rawp) |>
    pivot_longer(cols = c("c95", "c1", "rawp"), names_to = "identity", values_to = "num_seqs") |>
    mutate(identity = fct_relevel(identity, rev(c("c95", "c1", "rawp")))) |>
    mutate(
        cov = case_when(
            grepl("c3", sample) ~ "3X",
            grepl("c5", sample) ~ "5X",
            grepl("c10", sample) ~ "10X",
        ),
        dmg = case_when(
            grepl("mid", sample) ~ "mid",
            grepl("high", sample) ~ "high",
        ),
        assm = case_when(
            grepl("carpedeam-unsafe", assembler) ~ "Carpedeam\n(unsafe)",
            grepl("carpedeam-safe", assembler) ~ "Carpedeam\n(safe)",
            grepl("megahit", assembler) ~ "Megahit",
        )
    ) |>
    mutate(
        cov = fct_relevel(cov, c("3X", "5X", "10X")),
        dmg = fct_relevel(dmg, c("high", "mid"))
    )

# Plotting of Figure 3
colors <- c(
    "0.95" = "#3594cC",
    "1.0" = "#8cc5e3",
    "raw" = "#d8a6a6"
)

data |>
    mutate(identity = case_when(
        identity == "c95" ~ "0.95",
        identity == "c1" ~ "1.0",
        identity == "rawp" ~ "raw"
    )) |>
    mutate(identity = fct_relevel(identity, rev(c("0.95", "1.0", "raw")))) |>
    ggplot(aes(x = assm, y = num_seqs, fill = identity)) +
    geom_bar(stat = "identity", position = "stack", width = 1, color = "black", linewidth = 0.3) +
    scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3)) +
    labs(
        x = "Assembler",
        y = "# Proteins",
        fill = "Dereplication Level"
    ) +
    facet_wrap(~sample, scales = "free") +
    facet_grid2(c("dmg", "cov"),
        labeller = "label_both",
        scales = "free_y", independent = "y"
    ) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        legend.position = "bottom",
    ) +
    scale_fill_manual(values = colors, name = NULL)
# ) |>


# Now, let's check the degree of overlap between
files <- list.files("./data/genes/ancientGut", full.names = TRUE, pattern = "combined_cluster.tsv")

read_cluster_file <- function(filename) {
    pattern <- "^(.*)_(.*)-combined_cluster\\.tsv$"
    matches <- str_match(basename(filename), pattern)

    # Extract the matched groups
    sample <- matches[2]
    identity <- matches[3]

    df <- read_tsv(filename, col_names = c("rep", "member")) |>
        clean_names() |>
        mutate(sample = sample, identity = identity)
    return(df)
}

# Get the clustering of the protein of all datasets
clusters <- map_dfr(files, read_cluster_file)

clusters <- clusters %>%
    mutate(assembler = case_when(
        str_detect(member, "mh") ~ "mh",
        str_detect(member, "cp_s") ~ "cp_s",
        str_detect(member, "cp_u") ~ "cp_u",
        TRUE ~ "unknown"
    ))

# Get number of proteins
n_proteins <- clusters |>
    group_by(sample, identity) |>
    summarize(n = n(), .groups = "drop")

# Aggregate data at the cluster (rep) level
aggregated_data <- clusters %>%
    group_by(sample, identity, rep, assembler) %>%
    summarize(present = 1, .groups = "drop") %>%
    pivot_wider(names_from = assembler, values_from = present, values_fill = list(present = 0))



# Generate data suitable for an UpSet plot
# Just do one example
upset_data <- aggregated_data %>%
    filter(sample == "gut_sum_high_c10", identity == 0.95) %>%
    select(-rep, -sample, -identity) %>%
    as.data.frame()

# Plot the UpSet plot Figure 5
upset(upset_data, sets = colnames(upset_data), order.by = "freq", keep.order = TRUE)


# Find those unique for each assembler
unique_clusters <- aggregated_data |>
    group_by(sample, identity) |>
    filter(cp_s + cp_u + mh == 1) |>
    select(sample, identity, rep) |>
    inner_join(clusters, by = c("sample", "identity", "rep")) |>
    group_by(sample, identity, assembler) |>
    summarize(n = n(), .groups = "drop")

shared_all_clusters <- aggregated_data |>
    group_by(sample, identity) |>
    filter(cp_s + cp_u + mh == 3) |>
    select(sample, identity, rep) |>
    inner_join(clusters, by = c("sample", "identity", "rep")) |>
    group_by(sample, identity) |>
    summarize(n = n(), .groups = "drop")

carpedeam_clusters <- aggregated_data |>
    group_by(sample, identity) |>
    filter(cp_s == 1, cp_u == 1, mh == 0) |>
    select(sample, identity, rep) |>
    inner_join(clusters, by = c("sample", "identity", "rep")) |>
    group_by(sample, identity) |>
    summarize(n = n(), .groups = "drop")

other_clusters <- aggregated_data |>
    group_by(sample, identity) |>
    filter((cp_s == 1 & cp_u == 0 & mh == 1) | (cp_s == 0 & cp_u == 1 & mh == 1)) |>
    select(sample, identity, rep) |>
    inner_join(clusters, by = c("sample", "identity", "rep")) |>
    group_by(sample, identity) |>
    summarize(n = n(), .groups = "drop")

# Plot figure 5
data <- unique_clusters |>
    rename(category = assembler) |>
    bind_rows(shared_all_clusters |> mutate(category = "all")) |>
    bind_rows(carpedeam_clusters |> mutate(category = "cp")) |>
    bind_rows(other_clusters |> mutate(category = "other")) |>
    mutate(
        cov = case_when(
            grepl("c3", sample) ~ "3X",
            grepl("c5", sample) ~ "5X",
            grepl("c10", sample) ~ "10X",
        ),
        dmg = case_when(
            grepl("mid", sample) ~ "mid",
            grepl("high", sample) ~ "high",
        ),
        # category = case_when(
        #     grepl("cp-u", category) ~ "Carpedeam\n(unsafe)",
        #     grepl("cp-s", category) ~ "Carpedeam\n(safe)",
        #     grepl("mh", category) ~ "Megahit",
        #     TRUE ~ category
        # )
    ) |>
    mutate(
        cov = fct_relevel(cov, c("3X", "5X", "10X")),
        dmg = fct_relevel(dmg, c("high", "mid"))
    ) |>
    mutate(category = fct_relevel(category, rev(c("all", "cp", "other", "cp_s", "cp_u", "mh"))))



colors <- c(all = "#2066a8", cp = "#3594cC", other = "#8cc5e3", cp_s = "#a00000", cp_u = "#c46666", mh = "#d8a6a6")

ggplot(data, aes(x = identity, y = n, fill = category)) +
    geom_bar(stat = "identity", position = "stack", width = 0.9, color = "black", linewidth = 0.3) +
    scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3)) +
    labs(
        x = "Clustering identity",
        y = "#Proteins",
        fill = "Dereplication Level"
    ) +
    facet_wrap(~sample, scales = "free") +
    facet_grid2(c("dmg", "cov"),
        labeller = "label_both",
        scales = "free_y", independent = "y"
    ) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        legend.position = "bottom"
    ) +
    scale_fill_manual(values = colors, name = NULL) +
    guides(fill = guide_legend(nrow = 1))



# Let's see who this proteins map to
files <- list.files("./data/genes/ancientGut", full.names = TRUE, pattern = "_search.tsv")

read_search_file <- function(filename) {
    pattern <- "^(.*)_(.*)_(.*)_search\\.tsv$"
    matches <- str_match(basename(filename), pattern)

    # Extract the matched groups
    sample <- matches[2]
    assembler <- matches[3]
    identity <- matches[4]

    df <- read_tsv(filename, col_names = c("query", "target", "pident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits", "qlen", "tlen")) |>
        clean_names() |>
        mutate(sample = sample, assembler = assembler, identity = identity)
    return(df)
}

search_results <- map_dfr(files, read_search_file)

hits <- search_results |>
    mutate(qcov = alnlen / qlen) |>
    filter(pident >= 95) |>
    select(sample, assembler, identity, query) |>
    distinct()

# HOw many of the proteins have hits?
hits_proteins <- hits %>%
    group_by(sample, assembler, identity) %>%
    summarize(n = n(), .groups = "drop") |>
    inner_join(protein_stats |> select(sample, assembler, identity, num_seqs))


# Let's figure out if the proteins that are unique to each assembler have any hits
data <- aggregated_data |>
    group_by(sample, identity) |>
    filter(cp_s + cp_u + mh == 1) |>
    select(sample, identity, rep) |>
    inner_join(clusters, by = c("sample", "identity", "rep")) |>
    select(sample, identity, query = member) |>
    inner_join(hits, by = c("sample", "identity", "query")) |>
    group_by(sample, identity, assembler) |>
    summarize(hits = n(), .groups = "drop") |>
    inner_join(unique_clusters |>
        mutate(assembler = case_when(
            assembler == "cp_s" ~ "carpedeam-safe",
            assembler == "cp_u" ~ "carpedeam-unsafe",
            assembler == "mh" ~ "megahit"
        ))) |>
    mutate(no_hits = n - hits) |>
    select(sample, assembler, identity, hits, no_hits) |>
    pivot_longer(cols = c(hits, no_hits), names_to = "type", values_to = "n") |>
    mutate(
        cov = case_when(
            grepl("c3", sample) ~ "3X",
            grepl("c5", sample) ~ "5X",
            grepl("c10", sample) ~ "10X",
        ),
        dmg = case_when(
            grepl("mid", sample) ~ "mid",
            grepl("high", sample) ~ "high",
        ),
        assm = case_when(
            grepl("carpedeam-unsafe", assembler) ~ "Carpedeam\n(unsafe)",
            grepl("carpedeam-safe", assembler) ~ "Carpedeam\n(safe)",
            grepl("megahit", assembler) ~ "Megahit",
        )
    ) |>
    mutate(
        cov = fct_relevel(cov, c("3X", "5X", "10X")),
        dmg = fct_relevel(dmg, c("high", "mid")),
        type = fct_relevel(type, c("no_hits", "hits"))
    )

# Plot figures 6
colors <- c("hits" = "#424242", "no_hits" = "#AE4740")
ggplot(data |> filter(identity == "1.0"), aes(x = assm, y = n, fill = type)) +
    geom_bar(stat = "identity", position = "stack", width = 1, color = "black", linewidth = 0.3) +
    scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3)) +
    labs(
        x = "Assembler",
        y = "# Proteins",
    ) +
    facet_wrap(~sample, scales = "free") +
    facet_grid2(c("dmg", "cov"),
        labeller = "label_both",
        scales = "free_y", independent = "y"
    ) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        legend.position = "bottom",
    ) +
    scale_fill_manual(values = colors, name = NULL)

# Plot figures 7
ggplot(data |> filter(identity == "0.95"), aes(x = assm, y = n, fill = type)) +
    geom_bar(stat = "identity", position = "stack", width = 1, color = "black", linewidth = 0.3) +
    scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3)) +
    labs(
        x = "Assembler",
        y = "# Proteins",
    ) +
    facet_wrap(~sample, scales = "free") +
    facet_grid2(c("dmg", "cov"),
        labeller = "label_both",
        scales = "free_y", independent = "y"
    ) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        legend.position = "bottom",
    ) +
    scale_fill_manual(values = colors, name = NULL)


# Figure 8
data <- aggregated_data |>
    group_by(sample, identity) |>
    filter(cp_s + cp_u + mh == 3) |>
    select(sample, identity, rep) |>
    inner_join(clusters, by = c("sample", "identity", "rep")) |>
    select(sample, identity, query = member) |>
    left_join(hits, by = c("sample", "identity", "query")) |>
    mutate(assembler = ifelse(is.na(assembler), "no_hits", assembler)) |>
    group_by(sample, identity, assembler) |>
    summarize(hits = n(), .groups = "drop") |>
    pivot_wider(names_from = assembler, values_from = hits, values_fill = list(hits = 0)) |>
    mutate(hits = `carpedeam-safe` + `carpedeam-unsafe` + `megahit`) |>
    select(sample, identity, hits, no_hits) |>
    pivot_longer(cols = c(hits, no_hits), names_to = "type", values_to = "n") |>
    mutate(
        cov = case_when(
            grepl("c3", sample) ~ "3X",
            grepl("c5", sample) ~ "5X",
            grepl("c10", sample) ~ "10X",
        ),
        dmg = case_when(
            grepl("mid", sample) ~ "mid",
            grepl("high", sample) ~ "high",
        )
    ) |>
    mutate(
        cov = fct_relevel(cov, c("3X", "5X", "10X")),
        dmg = fct_relevel(dmg, c("high", "mid")),
        type = fct_relevel(type, c("no_hits", "hits"))
    )

colors <- c("hits" = "#424242", "no_hits" = "#AE4740")
ggplot(data, aes(x = identity, y = n, fill = type)) +
    geom_bar(stat = "identity", position = "stack", width = 1, color = "black", linewidth = 0.3) +
    scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3)) +
    labs(
        x = "Assembler",
        y = "# Proteins",
    ) +
    facet_wrap(~sample, scales = "free") +
    facet_grid2(c("dmg", "cov"),
        labeller = "label_both",
        scales = "free_y", independent = "y"
    ) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        legend.position = "bottom",
    ) +
    scale_fill_manual(values = colors, name = NULL)


# Figure 9
data <- aggregated_data |>
    group_by(sample, identity) |>
    filter(cp_s == 1, cp_u == 1, mh == 0) |>
    select(sample, identity, rep) |>
    inner_join(clusters, by = c("sample", "identity", "rep")) |>
    select(sample, identity, query = member) |>
    left_join(hits, by = c("sample", "identity", "query")) |>
    mutate(assembler = ifelse(is.na(assembler), "no_hits", assembler)) |>
    group_by(sample, identity, assembler) |>
    summarize(hits = n(), .groups = "drop") |>
    pivot_wider(names_from = assembler, values_from = hits, values_fill = list(hits = 0)) |>
    mutate(hits = `carpedeam-safe` + `carpedeam-unsafe`) |>
    select(sample, identity, hits, no_hits) |>
    pivot_longer(cols = c(hits, no_hits), names_to = "type", values_to = "n") |>
    mutate(
        cov = case_when(
            grepl("c3", sample) ~ "3X",
            grepl("c5", sample) ~ "5X",
            grepl("c10", sample) ~ "10X",
        ),
        dmg = case_when(
            grepl("mid", sample) ~ "mid",
            grepl("high", sample) ~ "high",
        )
    ) |>
    mutate(
        cov = fct_relevel(cov, c("3X", "5X", "10X")),
        dmg = fct_relevel(dmg, c("high", "mid")),
        type = fct_relevel(type, c("no_hits", "hits"))
    )

colors <- c("hits" = "#424242", "no_hits" = "#AE4740")
ggplot(data, aes(x = identity, y = n, fill = type)) +
    geom_bar(stat = "identity", position = "stack", width = 1, color = "black", linewidth = 0.3) +
    scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3)) +
    labs(
        x = "Assembler",
        y = "# Proteins",
    ) +
    facet_wrap(~sample, scales = "free") +
    facet_grid2(c("dmg", "cov"),
        labeller = "label_both",
        scales = "free_y", independent = "y"
    ) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        legend.position = "bottom",
    ) +
    scale_fill_manual(values = colors, name = NULL)


# Why if there are more contigs in CP we have the same number of proteins?

gff_files <- list.files("./data/genes/ancientGut", full.names = TRUE, pattern = "gff")

read_gff_files <- function(filename) {
    pattern <- "^(.*)_(.*)\\.gff$"
    matches <- str_match(basename(filename), pattern)

    # Extract the matched groups
    sample <- matches[2]
    assembler <- matches[3]

    df <- read_tsv(filename, col_names = FALSE, comment = "#") |>
        clean_names() |>
        select(contig_id = x1, start = x4, end = x5) |>
        mutate(sample = sample, assembler = assembler) |>
        mutate(ngene = row_number()) |>
        mutate(assm = case_when(
            assembler == "carpedeam-safe" ~ "cp_s",
            assembler == "carpedeam-unsafe" ~ "cp_u",
            assembler == "megahit" ~ "mh"
        )) |>
        mutate(prot_id = paste(sample, assm, ngene, sep = "-")) |>
        mutate(gene_length = end - start + 1) |>
        select(sample, assembler, assembler, contig_id, prot_id, gene_length)

    return(df)
}


# Get prodigal GFF files
gff_data <- map_dfr(gff_files, read_gff_files)

gff_data |>
    group_by(sample, assembler, contig_id) |>
    summarize(n = n(), coding_bp = sum(gene_length), .groups = "drop") |>
    inner_join(contigs) |>
    mutate(coding_density = coding_bp / length) |>
    filter(coding_density < 1 / 2) |>
    group_by(sample, assembler) |>
    summarise(n = n(), length = sum(length), .groups = "drop") |>
    mutate(
        cov = case_when(
            grepl("c3", sample) ~ "3X",
            grepl("c5", sample) ~ "5X",
            grepl("c10", sample) ~ "10X",
        ),
        dmg = case_when(
            grepl("mid", sample) ~ "mid",
            grepl("high", sample) ~ "high",
        ),
        assm = case_when(
            grepl("carpedeam-unsafe", assembler) ~ "Carpedeam\n(unsafe)",
            grepl("carpedeam-safe", assembler) ~ "Carpedeam\n(safe)",
            grepl("megahit", assembler) ~ "Megahit",
        )
    ) |>
    mutate(
        cov = fct_relevel(cov, c("3X", "5X", "10X")),
        dmg = fct_relevel(dmg, c("high", "mid")),
    ) |>
    ggplot(aes(x = n, y = length, fill = assm)) +
    geom_point(shape = 21, size = 3, color = "black") +
    # stat_summary(fun = median, geom = "line", size = 2, color = "#EB554A", group=1) +
    # stat_summary(fun = median, geom = "point", size = 2, fill = "#EB554A", shape = 21) +
    # scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3)) +
    labs(
        x = "# contigs",
        y = "Aggregated bp",
    ) +
    facet_grid2(c("dmg", "cov"),
        labeller = "label_both",
        scales = "free", independent = "all"
    ) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        legend.position = "bottom",
    ) +
    scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6)) +
    scale_fill_manual(values = c("#2D2D2D", "#EB554A", "#FFC300"), name = NULL)



# Get annotations for the contigs based on Prokka
prokka_gff_files <- list.files("./data/annotations/ancientGut", full.names = TRUE, pattern = "gff")

read_prokka_gff_files <- function(filename) {
    pattern <- "^(.*).raw-raw.proteins.(.*)\\.gff$"
    matches <- str_match(basename(filename), pattern)

    # Extract the matched groups
    sample <- matches[2]
    assembler <- matches[3]

    df <- read_tsv(filename, col_names = FALSE, comment = "#") |>
        clean_names() |>
        select(contig_id = x1, feature = x3, start = x4, end = x5) |>
        mutate(sample = sample, assembler = assembler) |>
        mutate(assembler = case_when(
            grepl("config7020", assembler) ~ "carpedeam-safe",
            grepl("config7022", assembler) ~ "carpedeam-unsafe",
            grepl("megahit", assembler) ~ "megahit"
        )) |>
        mutate(assm = case_when(
            assembler == "carpedeam-safe" ~ "cp_s",
            assembler == "carpedeam-unsafe" ~ "cp_u",
            assembler == "megahit" ~ "mh"
        )) |>
        mutate(feature_length = end - start + 1) |>
        select(sample, assembler, contig_id, feature, feature_length)

    return(df)
}

# Read the GFFs
prokka_annotations <- map_dfr(prokka_gff_files, read_prokka_gff_files)

prokka_features <- gff_data |>
    group_by(sample, assembler, contig_id) |>
    summarize(n = n(), coding_bp = sum(gene_length), .groups = "drop") |>
    inner_join(contigs) |>
    mutate(coding_density = coding_bp / length) |>
    # filter(coding_density < 1 / 2) |>
    inner_join(skani_search |> select(contig_id = query_name, sample, assembler, ani, align_fraction_query)) |>
    inner_join(prokka_annotations |> select(contig_id, sample, assembler, feature, feature_length))

# Which prokka features are in those contigs where less than 50% of their bps are not producing CDS
prokka_features |>
    filter(coding_density < 1 / 2) |>
    group_by(sample, assembler, feature) |>
    summarize(n = n(), length = sum(feature_length), .groups = "drop") |>
    # filter(feature != "CDS") |>
    mutate(
        cov = case_when(
            grepl("c3", sample) ~ "3X",
            grepl("c5", sample) ~ "5X",
            grepl("c10", sample) ~ "10X",
        ),
        dmg = case_when(
            grepl("mid", sample) ~ "mid",
            grepl("high", sample) ~ "high",
        ),
        assm = case_when(
            grepl("carpedeam-unsafe", assembler) ~ "Carpedeam\n(unsafe)",
            grepl("carpedeam-safe", assembler) ~ "Carpedeam\n(safe)",
            grepl("megahit", assembler) ~ "Megahit",
        )
    ) |>
    mutate(
        cov = fct_relevel(cov, c("3X", "5X", "10X")),
        dmg = fct_relevel(dmg, c("high", "mid")),
        feature = fct_reorder(feature, length)
    ) |>
    ggplot(aes(x = assm, y = length, fill = feature)) +
    geom_col(alpha = 0.8, color = "black") +
    labs(
        x = "Assembler",
        y = "Base pairs",
    ) +
    facet_grid2(c("dmg", "cov"),
        labeller = "label_both",
        scales = "free", independent = "all"
    ) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        legend.position = "bottom",
    ) +
    scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3)) +
    scale_fill_manual(values = c("#2D2D2D", "#EB554A", "#FFC300", "#F09868", "#91AEB7"), name = NULL)
