# from a cleaned list of variables defining the GALAXY_OUTPUT C struct
# generate numpy dtype with SHAPES preserved

BEGIN {
    print "# Automatically generated python structure to match the L-Galaxy binary output"
    print "import numpy"
    print "LGalaxiesStruct = numpy.dtype(["
    n = 0

    sizes["RNUM"] = RNUM
    sizes["NUM_METAL_CHANNELS"] = NUM_METAL_CHANNELS
    sizes["NUM_ELEMENTS"] = NUM_ELEMENTS
    sizes["NUM_COLDGAS_DUST_RATES"] = NUM_COLDGAS_DUST_RATES
    sizes["NUM_HOTGAS_DUST_RATES"] = NUM_HOTGAS_DUST_RATES
    sizes["SFH_NBIN"] = SFH_NBIN
}

function resolve_dim(x) {
    if (x in sizes) return sizes[x]
    return x + 0  # converts numeric strings safely, leaves symbols unchanged otherwise
}

{
    line = $0
    split(line, fields)

	if (fields[1] == "long" && fields[2] == "long") {
	    type = "long long"
	    name = fields[3]
	} else {
	    type = fields[1]
	    name = fields[2]
	}
	original_name = name

    # -----------------------------
    # dtype mapping
    # -----------------------------
    if (type == "float") {
        dtype = "numpy.float32"
        base_size = 4
    }
    else if (type == "long long") {
	    dtype = "numpy.int64"
	    base_size = 8
	}
    else if (type == "double" || type == "long") {
        dtype = "numpy.float64"
        base_size = 8
    }
    else if (type == "short") {
        dtype = "numpy.int16"
        base_size = 2
    }
    else {
        dtype = "numpy.int32"
        base_size = 4
    }

    # fix special cases
    if (type == "int") dtype = "numpy.int32"

    # -----------------------------
    # extract all dimensions
    # -----------------------------
    ndims = 0
    shape_str = ""

    n = split(line, tmp, /\[|\]/)

    for (i = 2; i <= n; i += 2) {
        dim = tmp[i]
        if (dim == "") continue

        dim = resolve_dim(dim)

        if (shape_str == "")
            shape_str = dim
        else
            shape_str = shape_str ", " dim

        ndims++
    }

    # -----------------------------
    # clean variable name
    # -----------------------------
    ia = match(name, /\[/)
    if (ia > 0)
        name = substr(name, 1, ia - 1)

    ia = index(name, ";")
    if (ia > 0)
        name = substr(name, 1, ia - 1)

    # -----------------------------
	# emit dtype entry (uniform tuple shapes)
	# -----------------------------
	
	if (ndims == 0) {
	    # scalar field → (1,)
	    print "('" name "', " dtype "),"
	}
	else if (ndims == 1) {
	    print "('" name "', " dtype ", (" shape_str ",)),"
	}
	else {
	    print "('" name "', " dtype ", (" shape_str ")),"
	}

    n++
}

END {
    print "('ending','i4',0)"
    print "])"

    print "properties_used = {}"
    print "for el in LGalaxiesStruct.names:"
    print "\tproperties_used[el] = True"
}