# check that data for an individual conforms to expected format
# by Karl Broman
check_individual <-
    function(ind, tol=1e-12)
{
    # list with two components, named "mat" and "pat"
    if(!is.list(ind) || length(ind)!=2 ||
       !all(names(ind) == c("mat", "pat")))
        stop('ind should be list with "mat" and "pat"')

    # check each chromosome
    for(i in 1:2) {
        # list with "alleles" and "locations" components
        if(!is.list(ind[[i]]) || length(ind[[i]])!=2 ||
           !all(names(ind[[i]]) == c("alleles", "locations")))
            stop('chromosome should be list with "alleles" and "locations"')

        alleles <- ind[[i]]$alleles
        locations <- ind[[i]]$locations

        # locations numeric
        if(!is.numeric(locations))
            stop("locations should be numeric")

        # alleles integer
        if(!is.integer(alleles))
            stop("alleles should be integers")

        # locations non-decreasing
        if( length(locations) > 1 && min(diff(locations)) < 0 )
            stop("locations should be non-decreasing")
    }

    # start and end positions are the same
    starts <- vapply(ind, function(a) a$location[1], 0.0)
    ends <- vapply(ind, function(a) a$location[length(a$location)], 0.0)
    stopifnot(abs(diff(starts)) < tol, abs(diff(ends)) < tol)

    TRUE
}
