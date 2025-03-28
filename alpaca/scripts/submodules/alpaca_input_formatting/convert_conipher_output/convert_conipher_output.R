library(jsonlite)
library(CONIPHER)
library(argparse)

parser = ArgumentParser(description = "Extract clone proportions and tree from CONIPHER output")
parser$add_argument("--CONIPHER_tree_object",type = "character",help = "path to CONIPHER tree object")
parser$add_argument("--output_dir",type = "character",help = "output directory")
args = parser$parse_args()


get_cp_table = function(alt_trees, alt_tree_id, clonality_table, ccf_cluster_table, trunk){
    cp_table = CONIPHER::compute_subclone_proportions(tree_list = alt_trees,
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