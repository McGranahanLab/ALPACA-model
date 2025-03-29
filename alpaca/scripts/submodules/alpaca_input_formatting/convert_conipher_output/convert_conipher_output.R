library(jsonlite)
library(optparse)

get_tree_level = function(tree_graph
                           , cluster)
{
  # trunk has level 1
  if (length(unique(as.numeric(tree_graph))) == 1){
    return(1)
  } else{
    trunk = unique(tree_graph[, 1][!tree_graph[, 1] %in% tree_graph[, 2]])
    if (cluster == trunk) return(1)
    else {
      clusters_in_tree = unique(as.numeric(as.matrix(tree_graph)))
      colnames(tree_graph) = c('Parent', 'Child')
      tree_graph = apply(tree_graph, c(1, 2), as.numeric)

      tree_df = as.data.frame(tree_graph)

      level = 1
      if (cluster %in% clusters_in_tree){
        current_cluster = as.numeric(cluster)
        while (current_cluster != trunk) {
          parent = tree_df[tree_df$Child == current_cluster, 'Parent']
          level = level + 1
          current_cluster = parent
        }
        return(level)
      } else return(NA)
    }
  }
}


compute_subclone_proportions = function(tree_list
                                         , ccf_cluster_table
                                         , clonality_table
                                         , trunk
                                         , force_clonal_100 = TRUE
                                         , tree_id = 1)
{
  options(stringsAsFactors = F)

  tree = tree_list[[tree_id]]

  clusters_in_tree = unique(c(as.matrix(tree)))

  # if tree is just clonal cluster
  if (length(clusters_in_tree) == 1){
    clone_proportion_table = (ccf_cluster_table > 0) * 100
  } else {

    colnames(tree) = c('Parent', 'Child')
    tree = as.data.frame(tree)
    region_IDs = colnames(ccf_cluster_table)

    # only use clusters in tree
    ccf_cluster_table = ccf_cluster_table[rownames(ccf_cluster_table) %in% clusters_in_tree, , drop = F]
    clonality_table = clonality_table[rownames(clonality_table) %in% clusters_in_tree, , drop = F]

    # force all clonal ccfs == 100 if flag is TRUE
    if (force_clonal_100 == T){
      clonal.factor = (clonality_table == 'clonal')
      ccf_cluster_table = ccf_cluster_table*(1 - clonal.factor) + 100*clonal.factor
    }

    # set a cap on CCFs above 100
    ccf_cluster_table[ccf_cluster_table > 100] = 100

    # assign clonal cluster to be == 100 ccf
    ccf_cluster_table[rownames(ccf_cluster_table) == trunk, ] = 100

    ccf_cluster_df = as.data.frame(ccf_cluster_table)
    ccf_cluster_df$cluster = rownames(ccf_cluster_df)


    ##### MAIN #####
    # Create a clone proportions dataframe to populate:
    proportions_df = as.data.frame(ccf_cluster_table)
    proportions_df$cluster = rownames(proportions_df)
    proportions_df[, region_IDs] = 0

    # Compute clone proportions for each region independently
    for (r in region_IDs){
      cols = c(r, 'cluster')
      region_ccf = ccf_cluster_df[, cols]
      colnames(region_ccf)[1] = c('ccf')

      clusters_present = region_ccf[region_ccf$ccf != 0, 'cluster']
      parents_present = unique(tree$Parent[tree$Parent %in% clusters_present])

      # order the parent subclones by tree level
      parent_df = data.frame(parent_node = parents_present)
      parent_df$level = sapply(parent_df$parent_node, function(p) return(get_tree_level(as.matrix(tree), p)))
      data.table::setorder(parent_df, level)

      # For each parent node in tree: compute the difference in CCF between parent node and sum of its children
      for (p in parent_df$parent_node){

        # TOP DOWN APPRAOCH:
        # Fix parental node CCFs, scale children nodes accordingly:
        # i.e. if children CCFs sum to > parental CCF, rescale child CCFs to be proportional values of parent ccf (so that their CCFs sum to parent exactly)

        children_nodes = tree[tree$Parent == p, 'Child']
        parent_ccf = as.numeric(region_ccf[region_ccf$cluster == p, 'ccf'])
        sum_children_ccf = sum(region_ccf[region_ccf$cluster %in% children_nodes, 'ccf'])

        if (sum_children_ccf > parent_ccf){
          parent_proportion = 0
          region_ccf[region_ccf$cluster %in% children_nodes, 'ccf'] = parent_ccf * region_ccf[region_ccf$cluster %in% children_nodes, 'ccf'] / sum(region_ccf[region_ccf$cluster %in% children_nodes, 'ccf'])
        } else {
          parent_proportion = parent_ccf - sum_children_ccf
        }
        proportions_df[proportions_df$cluster == p, r] = parent_proportion

        # if the clones are terminal, add them to proportions_df as well
        for (d in children_nodes){
          if (d %in% tree[, 'Child'] & !(d %in% tree[, 'Parent'])) {
            proportions_df[proportions_df$cluster == d,r] = region_ccf[region_ccf$cluster == d, 'ccf']
          }
        }
      }
    }

    # fix final output matrix
    rownames(proportions_df) = proportions_df$cluster
    clone_proportion_table = as.matrix(proportions_df[, region_IDs, drop = F])
  }
  return(clone_proportion_table)
}

get_cp_table = function(alt_trees, alt_tree_id, clonality_table, ccf_cluster_table, trunk){
    cp_table = compute_subclone_proportions(tree_list = alt_trees,
                                                    ccf_cluster_table = ccf_cluster_table,
                                                    clonality_table = clonality_table,
                                                    trunk = trunk,
                                                    force_clonal_100 = TRUE,
                                                    tree_id = alt_tree_id)
    cp_table = cp_table / 100
    cp_table = as.data.frame(cp_table)
    cp_table$clone = paste0('clone', rownames(cp_table))
    colnames(cp_table) = gsub('\\.', '-', colnames(cp_table))
    return (cp_table)
}


extractTreeGraphPaths = function(tree_graph){
    # function that extracts all end-to-end tree paths for a tree graph
    clones_in_tree = unique(as.numeric(as.matrix(tree_graph)))
    if (length(clones_in_tree)==1) paths = list(1)
    else {
        tree_graph = apply(tree_graph, c(1,2), as.numeric)
        tree_graph = as.data.frame(tree_graph)
        colnames(tree_graph) = c('Parent', 'Child')
        paths = list()
        trunk = unique(tree_graph$Parent[!tree_graph$Parent %in% tree_graph$Child])
        terminal.clones = tree_graph[!(tree_graph$Child %in% tree_graph$Parent), 'Child']
        for (terminal in terminal.clones){
        p = c(terminal)
        current_clone = terminal
        while (current_clone != trunk){
            parent = tree_graph[tree_graph$Child==current_clone, 'Parent']
            p = c(p,parent)
            current_clone = parent  
        }
        p = rev(p)
        paths = append(paths, list(p))
        }
        return(paths)
    }
  }



# Create option list
option_list = list(
  make_option(c("--CONIPHER_tree_object"), type="character", 
              help="path to CONIPHER tree object"),
  make_option(c("--output_dir"), type="character", 
              help="output directory")
)

# Create the parser
parser = OptionParser(
  description = "Extract clone proportions and tree from CONIPHER output",
  option_list = option_list
)

# Parse the arguments
args = parse_args(parser)

# Read the tree object and set output directory
tree_object = readRDS(args$CONIPHER_tree_object)
output_dir = args$output_dir

default_tree = tree_object$graph_pyclone$default_tree
tree_paths = extractTreeGraphPaths(tree_graph = default_tree)
tree_path_clone_names = lapply(tree_paths, function(path){paste0('clone', path)})
tree_json_path = sprintf('%s/tree_paths.json',output_dir)
write(toJSON(tree_path_clone_names), tree_json_path)

alt_trees = tree_object$graph_pyclone$alt_trees
clonality_table = tree_object$clonality_out$clonality_table_corrected
ccf_cluster_table = tree_object$nested_pyclone$ccf_cluster_table
trunk = tree_object$graph_pyclone$trunk 
cp_table = get_cp_table(alt_trees, 1, clonality_table, ccf_cluster_table, trunk)
cp_table_path = sprintf('%s/cp_table.csv',output_dir)
write.csv(cp_table, cp_table_path, row.names = FALSE)
print('Done')