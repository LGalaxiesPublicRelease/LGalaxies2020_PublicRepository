# extract a block starting with a line containing GALAXY_OUTPUT and ending with a }
BEGIN {
flag=0
line=""
d=""
ia=0
b=""
}
/GALAXY_OUTPUT/,/\}/ {

    line=$0

    if(flag == 0 && line ~ /\{/) {
        flag=1
        next
    }

    if(line ~ /\}/)
        exit

    if(flag != 1)
        next

    # Skip preprocessor directives
    if(line ~ /^[ \t]*#/)
        next

    # Skip block comments
    if(line ~ /^[ \t]*\/\*/)
        next

    if(line ~ /^[ \t]*\*/)
        next

    # Skip blank lines
    if(line ~ /^[ \t]*$/)
        next

    # Remove trailing // comment
    sub(/[ \t]*\/\/.*/, "", line)

    # Trim whitespace
    gsub(/^[ \t]+|[ \t]+$/, "", line)

    # Remove trailing semicolon
    sub(/;$/, "", line)

    if(line != "")
        print line
}