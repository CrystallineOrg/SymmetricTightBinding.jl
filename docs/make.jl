using Documenter
using SymmetricTightBinding
using Optim # load the Optim extension so its docstrings (`fit`, …) are available to `@docs`

# use SymmetricTightBinding as the default module for doctests, and also load Crystalline
DocMeta.setdocmeta!(
    SymmetricTightBinding,
    :DocTestSetup,
    :(using SymmetricTightBinding, Crystalline);
    recursive = true,
)

makedocs(;
    sitename = "SymmetricTightBinding.jl",
    authors = "Antonio Morales Perez <antonio.morales@dipc.org>, Thomas Christensen <thomas@dtu.dk>, and contributors",
    # `fit` & `multistart_fit` are stubs in the main package, with their methods and
    # docstrings attached in the Optim extension; that module must be listed here for
    # Documenter to find those docstrings
    modules = [
        SymmetricTightBinding,
        Base.get_extension(SymmetricTightBinding, :SymmetricTightBindingOptimExt),
    ],
    repo = Remotes.GitHub("CrystallineOrg", "SymmetricTightBinding.jl"),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://CrystallineOrg.github.io/SymmetricTightBinding.jl",
        size_threshold = 1000000,
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Band symmetry" => "band-symmetry.md",
        "Berry curvature & Chern numbers" => "berry.md",
        "Density of states" => "dos.md",
        "Symmetry breaking" => "symmetry-breaking.md",
        "Non-Hermitian models" => "nonhermitian.md",
        "Fitting to band structures" => "fitting.md",
        "Designing models to an objective" => "design.md",
        "API" => "api.md",
        "Internal API" => "internal-api.md",
        "Theory" => "theory.md",
        "Developer notes" => [
            "Overview" => "devdocs/README.md",
            "devdocs/trs_notes.md",
            "devdocs/fourier.md",
            "devdocs/1d_example.md",
            "devdocs/symmetry_eigenvalue_conventions.md",
        ],
    ],
    warnonly = Documenter.except(
        :autodocs_block,
        :cross_references,
        :docs_block,
        :doctest,
        :eval_block,
        :example_block,
        :footnote,
        :linkcheck_remotes,
        :linkcheck,
        :meta_block,
        :parse_error,
        :setup_block,
        #:missing_docs # don't fail on missing doc strings, too annoying
    ),
    clean = true,
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(;
    repo = "github.com/CrystallineOrg/SymmetricTightBinding.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
    push_preview = true, # deploy docs for PRs
)