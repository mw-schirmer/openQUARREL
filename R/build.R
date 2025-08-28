setwd("D:/mygithub/openQUARREL")
pkgdown::build_site()
usethis::use_pkgdown_github_pages()
pkgdown::deploy_to_branch()
