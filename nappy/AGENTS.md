# Guidlines for nappy development

## Python coding rules

- Write code targeting Python 3.9 or later.
- Use the docopt package for standalone scripts that take command arguments and options.
- Make `main()` function instead of writing code just under `if __name__ == "__main__":`.
- Use `python_script_template.py` as a template for script. If it is not used as a CLI tool, you don't need to use this template.
- Set YYMMDD to `__version__`, where YYMMDD is the creation or update date.

## matplotlib

- Use DPI greater than or equal 150.
- Make the aspect ratio 1:1 if possible or appropriate, regardless the physical units of x and y axes.
- Use seaborn to make the style consistent as,
  ```python
  import seaborn as sns
  sns.set_theme(context='talk', style='ticks')
  ```
