# Guidelines for AI Agents

Rules for AI agents modifying the nap package.
The nap package contains a molecular dynamics program called pmd and a Python pre/post-processing package called nappy.
The following directories at the root level may be ignored:
- `JOSS_paper`
- `mkconf`
- `neb`
- `not_used`
- `qmcl`


## Programming Language Rules

### pmd

- Source code for the Fortran program pmd goes under `pmd/`.
- Fortran code should conform to Fortran90. Fortran 2002 features may be used when necessary.
- Indent with 2 spaces. Do not use tabs.
- Keep lines to 78 characters where possible; wrap long lines with `&` at the end.

### nappy

- The Python package nappy lives under `nappy/`.
- Write code targeting Python 3.9 or later.


## Git Rules

- Work on the `ai-dev` branch.
- Obtain permission before running `git commit`.
- Commit in logical units (one feature addition or one bug fix per commit). Do not batch many unrelated changes into a single commit.
- Use Conventional Commits for commit messages.

## Work Reporting

- When implementing or modifying anything, save the work plan to `agents/dev/plan_YYMMDD_NAME.md`, the todo list to `agents/dev/todo_YYYMMDD_NAME.md`, and the work log to `agents/dev/log_YYMMDD_NAME.md`, where NAME is a brief description of the work and YYMMDD is the date work began.
