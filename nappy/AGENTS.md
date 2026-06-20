# Guidlines for nappy development

## Python coding rules

- Write code targeting Python 3.9 or later.
- Use the docopt package for standalone scripts that take command arguments and options.
- Make `main()` function instead of writing code just under `if __name__ == "__main__":`.
- Use `python_script_template.py` as a template for script. If it is not used as a CLI tool, you don't need to use this template.
- Set YYMMDD to `__version__`, where YYMMDD is the creation or update date.
