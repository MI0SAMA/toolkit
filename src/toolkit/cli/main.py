import click
from toolkit.cli.convert import convert
from toolkit.cli.analyze import analyze
from toolkit.cli.manipulate import manipulate


@click.group()
@click.version_option(version='0.1.0')
def main():
    """Toolkit - Theoretical calculation toolkit.

    Run 'tk' without arguments for interactive menu.
    """
    pass


main.add_command(convert)
main.add_command(analyze)
main.add_command(manipulate)


def _entry():
    """CLI entry point for interactive menu fallback."""
    import sys
    if len(sys.argv) == 1:
        from toolkit.cli.menu import run_menu
        run_menu()
    else:
        main()


if __name__ == '__main__':
    _entry()
