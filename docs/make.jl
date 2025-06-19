using Documenter
using SymmetricTightBinding

makedocs(
    sitename = "SymmetricTightBinding.jl",
    authors = "Antonio Morales Perez <antonio.morales@dipc.org>, Thomas Christensen <thomas@dtu.dk>, and contributors",
    modules = [SymmetricTightBinding],
    repo = Remotes.GitHub("CrystallineOrg", "SymmetricTightBinding.jl"),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://CrystallineOrg.github.io/SymmetricTightBinding.jl",
        size_threshold = 1000000
    ),
    pages = [
        "Home"                  => "index.md",
        "API"                   => "api.md",
        "Internal API"          => "internal-api.md",
    ],
    warnonly = Documenter.except(
        :autodocs_block, :cross_references, :docs_block, :doctest, :eval_block, 
        :example_block, :footnote, :linkcheck_remotes, :linkcheck, :meta_block, 
        :parse_error, :setup_block,
        #:missing_docs # don't fail on missing doc strings, too annoying
    ),
    clean = true
)

# use SymmetricTightBinding as the default module for doctests
DocMeta.setdocmeta!(
    SymmetricTightBinding,
    :DocTestSetup,
    :(using SymmetricTightBinding);
    recursive=true
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo   = "github.com/CrystallineOrg/SymmetricTightBinding.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
    push_preview = true # deploy docs for PRs
)