import pep8
fchecker = pep8.Checker('snafu.py', show_source=True)
file_errors = fchecker.check_all()

