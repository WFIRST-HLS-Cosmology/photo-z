# print output to standard error stream

errlog() { 
    printf "$@\n" 1>&2
}
