######
# title: dependency.py
#
# language: python3
# date: 2026-08-14
# license: BSD-3-Clause
# author: Elmar Bucher, Heber Rocha, Claude Anthropic
#
# description:
#   pcdl, by default, only installs and imports lightweight core dependencies.
#   the code below deals with installation and import of the more specialized
#   heavyweight libraries.
######


# library
import importlib


# function
def optional_import(s_module, s_attr=None, s_pip=None, s_caller=None):
    """
    input:
        s_module: string
            dotted module path to import, e.g. 'anndata' or 'bioio.writers'.

        s_attr: string; default None
            if given, equivalent to from s_module import s_attr.

        s_pip: string; default None
            pip install name for the module, if this differs from s_module
            (e.g. s_module='skimage' but s_pip='scikit-image').
            if None, the first dot-separated part of s_module is used.

        s_caller: string; default None
            name of the calling pcdl function, to mention in the error message.

    output:
        the imported module, or, if s_attr is given, the requested attribute
        of the imported module.

    description:
        function lazily load an optional pcdl dependency.
        if the library is not installed, an error message is raised,
        pointing the user to install the missing library.
    """
    s_pip = s_module.split('.')[0] if (s_pip is None) else s_pip
    try:
        o_module = importlib.import_module(s_module)
    except ImportError as e:
        s_fct = s_caller if not (s_caller is None) else ' '
        raise ModuleNotFoundError(
            f"Error{s_fct}: this functionality requires the optional dependency '{s_pip}', which is not installed.\n" +
            f"pcdl was installed light weight (default), without this and other heavyweight, specialized libraries.\n" +
            f"to fix this, either:\n" +
            f"+ install the missing library manually: pip install {s_pip}\n" +
            f"+ or install pcdl with all optional dependencies: pip install pcdl[full]\n"
        ) from e
    return o_module if (s_attr is None) else getattr(o_module, s_attr)
