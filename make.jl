push!(LOAD_PATH,"src/")

using Pkg;

Pkg.instantiate()

using Documenter, HighsDocs

# makedocs(sitename="HiGHS Documentation",format = Documenter.HTML(
# ))


makedocs(
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        # prettyurls = !("local" in ARGS),
        prettyurls = get(ENV, "CI", nothing) == "true",
        highlights = ["yaml"],
    ),
    clean = false,
    sitename = "HiGHS Documentation",
    authors = "Julian Hall and Ivet Galabova",
    pages = [
        "About" => "index.md",
        "Guide" => "guide.md", 
        "HiGHS in Python" => Any[
            "Get started in Python" => "python/pip.md",
            "Enums" => "python/enums.md",
            "Classes" => Any[
	        "Introduction" => "python/classes/Index.md",
	        "HighsSparseMatrix" => "python/classes/HighsSparseMatrix.md",
	        "HighsLp" => "python/classes/HighsLp.md",
	        "HighsSolution" => "python/classes/HighsSolution.md",
	        "HighsBasis" => "python/classes/HighsBasis.md",
	        "HighsInfo" => "python/classes/HighsInfo.md",
	        "Other" => "python/classes/Other.md"],
            "Examples" => "python/example-py.md",
            "Notebooks" => "python/notebooks.md"],
        "HiGHS in C++" => Any[
            "Get started in C++" => "cpp/get-started.md",
            "The HiGHS library" => "cpp/library.md",
            "Linking" => "cpp/link.md",
            "Examples" => "cpp/examples.md"],
        "Binaries" => "binaries.md",
        "Executable" => "executable.md",
	"Options" => Any[
	    "Introduction" => "options/intro.md",
	    "Definitions" => "options/definitions.md"],
        "Interfaces" => "interfaces.md",
        "Terminology" => "terminology.md",
],
    strict = !("strict=false" in ARGS),
    doctest = ("doctest=only" in ARGS) ? :only : true,
)