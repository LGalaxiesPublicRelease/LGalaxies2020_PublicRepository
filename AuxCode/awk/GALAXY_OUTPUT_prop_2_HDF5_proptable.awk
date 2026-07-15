# from a cleaned list of variables defining the GALAXY_OUPUT C struct, create a .h file containing the
# LGalaxy struct required for HDF5 output.
BEGIN {
    # recognised C declaration types
    types["float"]=1
    types["double"]=1
    types["int"]=1
    types["short"]=1
    types["long"]=1
    types["char"]=1
    types["struct"]=1
}

{
    line = $0

    # trim leading/trailing whitespace
    gsub(/^[ \t]+|[ \t]+$/, "", line)

    # ignore blank lines
    if (line == "")
        next

    # ignore comments
    if (line ~ /^\/\//)
        next
    if (line ~ /^\/\*/)
        next
    if (line ~ /^\*/)
        next

    split(line, fields)

    type = fields[1]

    # only accept genuine C declarations
    if (!(type in types))
        next

    if (type == "struct")
        fieldname = fields[3]
    else if (type == "long" && fields[2] == "long")
        fieldname = fields[3]
    else
        fieldname = fields[2]

    # remove array dimensions
    sub(/\[.*/, "", fieldname)

    # remove trailing semicolon
    sub(/;.*/, "", fieldname)

    unit = ""
    desc = ""

    if (match(line, /\/\/[^\/]*/)) {
        unit = substr(line, RSTART+2, RLENGTH-2)

        rest = substr(line, RSTART+RLENGTH)

        if (match(rest, /\/\/.*/))
            desc = substr(rest, RSTART+2)
    }

    gsub(/^[ \t]+|[ \t]+$/, "", unit)
    gsub(/^[ \t]+|[ \t]+$/, "", desc)

    print fieldname " , " unit " , " desc
}