extract_link = function(fit, par = 'gamma'){
    model = get_stancode(fit)
    stringr::str_split(regmatches(model, regexpr(sprintf('//%s_link[a-z_:]+', par), model)), ':', simplify = T)[1,2]
}
