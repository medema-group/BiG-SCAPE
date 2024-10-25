""" Configures cli """

import click


class CommandOrder(click.Group):
    """Modification of click Group to conserve command order"""

    def list_commands(self, ctx: click.Context) -> list[str]:
        """list all sub-commands in the order they were supplied

        Args:
            ctx (click.Context): Click context

        Returns:
            list[str]: names of sub-commands
        """
        return list(self.commands)
