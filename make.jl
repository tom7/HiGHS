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
    authors = "Julian Hall, Ivet Galabova and Michael Feldmeier.",
    pages = [
        "Home" => "index.md",
        "About" => "about.md", 
        "Guide" => "guide.md", 
        "HiGHS in Python" => Any[
            "Get Started in Python" => "python/pip.md",
            "Example" => "python/example-py.md",
            "Notebooks" => "python/notebooks.md"],
        "HiGHS in C++" => Any[
            "Get Started in C++" => "cpp/get-started.md",
            "The HiGHS Library" => "cpp/library.md",
            "Options" => "cpp/options.md",
            "Linking" => "cpp/link.md",
            "Examples" => "cpp/examples.md"],
        "Interfaces" => "interfaces.md",
        "Terminology" => "terminology.md",
],
    strict = !("strict=false" in ARGS),
    doctest = ("doctest=only" in ARGS) ? :only : true,
)