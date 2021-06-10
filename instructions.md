# Assignments

Edition 2020-2021

This file gives a detailed overview of what you have to do for this project. 

## In brief

For the exam project, you pick a optimization related topic of your interest that is *not* covered in detail in class. This can be an algorithm, an application you solve with methods seen the course or some theoretical aspect you want to study. You write some code that you add to the STMOZOO codebase (including documentation, tests etc.) and illustrate you application in a notebook.

## Getting started

- [x] pick a project (take a look at `project ideas.md` or discuss with Michiel)
- [x] [fork](https://docs.github.com/en/enterprise-server@2.20/github/getting-started-with-github/fork-a-repo) this repo
- [x] create a new branch with a short indicative name, e.g., `GeneticProgramming`. **Don't use spaces in the name!**
- [x] make a local clone of the repository 
- [x] open a [pull request](https://docs.github.com/en/free-pro-team@latest/desktop/contributing-and-collaborating-using-github-desktop/creating-an-issue-or-pull-request) to the **master** branch of this repo. This makes it clear you are starting the project!

## Source code

Every project needs to have some source code, at least one function! You have to decide which parts belong in the source code (and can hence be readily loaded by other users) and which parts of your project will be in the notebook where people can see and interact with your code.

Developing code can be done in any text editor, though we highly recommend [Visual Studio Code](https://code.visualstudio.com/), with Juno the environment for Julia. [Atom](https://atom.io/) is an alternative but is not supported anymore. When developing, you have to activate your project. Assuming that the location of the REPL is the project folder, open the Pkg manager (typing `]`) and type `activate .`. The dot indicated the current directory. If you use external packages in your project, for example, Zygote or LinearAlgebra, you have to add them using `add PACKAGE` in the package manager. This action will create a dependency and update the `Project.toml` file.

Importantly, all your code should be in a [module](https://docs.julialang.org/en/v1/manual/modules/), where you export only the functions useful for the user.

- [x] In the `src` folder, add a new Julia file with your source code, for example `geneticprogramming.jl`. Don't use spaces or capitals in the file name.
- [x] Link your file in `STMOZOO.jl` using `include(filename)`,  running the code.
- [x] Create a module environment in your file for all your code. Use [camel case](https://en.wikipedia.org/wiki/Camel_case) for the name.
  - use `module GeneticProgramming begin ... end` to wrap your code;
  - import everything you need from external packages: `using LinearAlgebra: norm`;
  - export your functions using `export`
- [x] write awesome code!
- [x] take a look at your code regarding the [Julia style guide](https://docs.julialang.org/en/v1/manual/style-guide/)
- [x] check the [Julia performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/)
- [x] document *every* function! Make sure that an external user can understand everything! Be liberal with comments in your code. Take a look at the [guidelines](https://docs.julialang.org/en/v1/manual/documentation/)

## Unit tests

Great, we have written some code. The question is, does it work? Likely you have experimented in the REPL. For a larger project, we would like to have guarantees that it works, though. Luckily, this is very easy in Julia, where we can readily use [Unit testing](https://docs.julialang.org/en/v1/stdlib/Test/).

You will have to write a file with some unit tests, ideally testing every function you have written! The fraction of functions that are tested is called [code coverage](https://en.wikipedia.org/wiki/Code_coverage). This project is monitored automatically using Travis (check the button on the readme page!). Currently, coverage is 100%, so help to keep this as high as possible!

Tests can be executed using the `@test` macro. You evaluate some functions and check their results. The result should evaluate to `true`. For example: `@test 1+1 == 2` or `@test √(9) ≈ 3.0`. 

It makes sense to group several tests, which can be done using `@testset "names of tests" begin ... end`.

Your assignments:
- [x] add a source file to the `test/` folder, the same name as your source code file.
- [x] add an `include(...)` with the filename in `runtests.jl`
- [x] in your file, add a block `@testset "MyModule" begin ... end` with a series of sensible unit tests. Use subblocks of `@testset` if needed.
- [x] run your tests, in the package manager, type `test`. It will run all tests and generate a report.

Travis will automatically run your unit tests online when you push to the origin repo.

## Documentation

Hopefully, you have already documented all your functions, so this should be a breeze! We will generate a documentation page using the [Documenter](https://juliadocs.github.io/Documenter.jl/stable/man/guide/) package. Since we will not put the project in the package manager, we won't host the documentation, though we generate HTML pages anyway.

- [x] add markdown file to `docs/src/man` with the documentation.
- [x] write a general introduction explaining the rationale of your code.
- [x] use a `@docs` block to add your functions with their documentation.
- [x] update the `make.jl` file, linking your page.
- [x] run the `make.jl` file to generate the documentation, an HTML file, not added to the repo.

## Notebook

Finally, you have to add a [Pluto](https://github.com/fonsp/Pluto.jl) notebook to the `notebook` folder. Again use the same name you used for your source code. Depending on the nature of your project, this will be the most extensive task! Make full use of Pluto's interactivity to illustrate your code. In contrast to the documentation page, this is not the place to explain your functions but rather show what you can do with your software or explain a concept.

## Code review

Each of you will have to perform a code review of two other projects on **January 7**. You have the full day to do this, though it should not take too long. The aim is to **help** the other groups to make each other's project even better.

- [ ] clone the projects you have to review, so you can run them locally;
- [ ] check the source code, is the documentation clear? Anything obvious that can be improved.
- [ ] run the tests. Do they work? Anything that could be tested but is not done so?
- [ ] Is the documentation page clear? Do you find any typos? Could an example be added?
- [ ] Take a look at the notebook. Any suggestions there to improve this?

You can leave comments in the pull request page. If you have things you want to fix yourself you can also perform a pull request like I did in your code review. **Communicate to the person when you are finised.**

on **January 8**, you have the full day to:
- [ ] merge the entire request and fix any issues you find meaningful.
- [ ] fill in a small questionnaire on Ufora about your project and the projects you have reviewed.

Please address all the comments you get as a common courtesy. This does not mean you have to implement all suggestions, but you do have to clarify on why you won't/can't do them.

Afterward, I will approve all pull requests, finalizing the package.

> Two days are blocked for the code review because some of you are in different time zones. This does **not** mean you will need these days in full. You will need approximately reserve one hour for each day. For the remainder, you can do whatever you want.
