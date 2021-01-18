# Rhub test script
require(rhub)
## validate_email()

## start cran like testing, save for report
cran_prep <- check_for_cran()

### add report to cran-comments.md
write(cran_prep$cran_summary(), file="cran-comments.md", append = FALSE, sep = " ")

## evaluate results
### get checks

previous_checks <- rhub::list_package_checks("~/Gits/consseg",
                                             email = "fall@bioinf.uni-leipzig.de",
                                             howmany = 4)

### list checks
previous_checks

### show checks
#### for i in 0:len(previous_checks$group)
group_id <- previous_checks$group[1]
group_check <- rhub::get_check(group_id)
group_check

