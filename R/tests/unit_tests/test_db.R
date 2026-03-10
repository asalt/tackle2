library(here)
library(testthat)
library(DBI)
library(RSQLite)
library(stringr)

TESTDB <- tempfile("tackle2-test-", fileext = ".db")

sim_tools <- new.env()
source(file.path(here("R"), "simulate.R"), local = sim_tools)

io_tools <- new.env()
source(file.path(here("R"), "./io.R"), local = io_tools)

db_tools <- new.env()
source(file.path(here("R"), "db.R"), local = db_tools)

fgsea_tools <- new.env()
source(file.path(here("R"), "fgsea.R"), local = fgsea_tools)

geneset_tools <- new.env()
source(file.path(here("R"), "geneset_utils.R"), local = geneset_tools)


testthat::test_that("test db table setup", {
  print(TESTDB)
  on.exit(if (file.exists(TESTDB)) unlink(TESTDB), add = TRUE)
  db_tools$initialize_db(db_path=TESTDB)
  con <- db_tools$get_con(TESTDB)
  on.exit(db_tools$close_con(con), add = TRUE)
  res <- dbGetQuery(con, "SELECT * FROM sqlite_master WHERE type='table';")
  for (table in c("Ranks", "RankObjects", "Collections", "Pathways", "CollectionResults", "Curves")){
    expect_true(table %in% res$name, info = table)
  }
}
)


setup <- function() {

  print(TESTDB)
  db_tools$initialize_db(db_path=TESTDB)
  con <- db_tools$get_con(TESTDB)


  data1 <- sim_tools$simulate_preranked_data(seed = 1234, sample_frac = .4)
  data2 <- sim_tools$simulate_preranked_data(seed = 4321, sample_frac = .4)
  data <- list(first = data1, second = data2)

  test_data <- sim_tools$generate_test_data(pathways = c("H", "GO:BP"), preranked_data = data)

  genesets_info <- list(
                    "H_" = geneset_tools$get_collection(category="H", subcategory = ""),
                    "C5_GO:BP" = geneset_tools$get_collection(category="C5", subcategory = "GO:BP")
                    )
  genesets_list_of_lists <- purrr::map(genesets_info, geneset_tools$genesets_df_to_list)

  collections <- names(test_data)
  rank_names <- names(test_data[[1]])

  # INSERT rank names
  #for (rank_name in rank_names){
  for (ix in seq_along(rank_names)){
      rankobj_id <- db_tools$insert_rankobj(con = con, rank_name = rank_names[[ix]])
      db_tools$insert_ranks(con = con, rankobj_id = rankobj_id, ranks_data = data[[ix]])
  }

  # INSERT collections
  for (collection in collections){ # this is slow and could be optimized
      collection_id <- db_tools$insert_collection(con, collection)
      for (pathway_name in names(genesets_list_of_lists[[collection]])) {
        db_tools$insert_pathway(
            con,
            collection_id=collection_id,
            pathway_name = pathway_name,
            members = str_c(genesets_list_of_lists[[collection]][[pathway_name]])
                )
       }
  }


  # INSERT gsea results
  for (rank_name in rank_names){
    for (collection in collections){
      db_tools$insert_results(con, rank_name = rank_name, collection_name = collection,
      results = test_data[[collection]][[rank_name]]
      )
    }
  }


  rankobjs <- io_tools$ranks_dfs_to_lists(data)

  db_tools$insert_ranks(con, rankobj_name = names(rankobjs)[[1]], ranks_data = rankobjs[[1]])
  db_tools$insert_ranks(con, rankobj_name = names(rankobjs)[[2]], ranks_data = rankobjs[[2]])

  rankorder1 <- fgsea_tools$get_rankorder(
    genesets_list_of_lists[[1]][[1]],
    rankobjs[[1]]
  )

  rankorder2 <- fgsea_tools$get_rankorder(
    genesets_list_of_lists[[1]][[1]],
    rankobjs[[2]]
  )



  # INSERT
  db_tools$insert_curve(
    con,
    rankobj_name = names(rankobjs)[[1]],
    pathway_name = names(genesets_list_of_lists[[1]])[[1]],
    curve_data = rankorder1$curve
  )

  # INSERT
  db_tools$insert_curve(
    con,
    rankobj_name = names(rankobjs)[[1]],
    pathway_name = names(genesets_list_of_lists[[1]])[[2]],
    curve_data = rankorder1$curve
  )

  # INSERT
  db_tools$insert_curve(
    con,
    rankobj_name = names(rankobjs)[[2]],
    pathway_name = names(genesets_list_of_lists[[1]])[[1]],
    curve_data = rankorder1$curve
  )

  # INSERT
  db_tools$insert_curve(
    con,
    rankobj_name = names(rankobjs)[[2]],
    pathway_name = names(genesets_list_of_lists[[1]])[[2]],
    curve_data = rankorder1$curve
  )





    # pathway_name = names(genesets_list_of_lists[[1]])[[1]]

}
