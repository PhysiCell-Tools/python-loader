#os.system(f'python3 -c "import pcdl; help(pcdl.TimeStep.get_multicellds_version)" >> {s_opathfile}')


# library
import pcdl
import os


# function
def docstring_md(s_function, ls_doc, s_header=None, s_opath='./docstring/'):
    os.makedirs(s_opath, exist_ok=True)
    s_opathfile = s_opath + f'{s_function}.md'
    f = open(s_opathfile, 'w')
    if (s_header is None):
        s_header = '{s_function}()'
    f.write('# {s_header}\n\n')
    for s_doc in ls_doc:
        if (len(s_doc.strip()) == 0):
            pass
        elif (s_doc.find('input:') > -1): 
            f.write(f'## {s_doc.strip()}\n```\n')
        elif (s_doc.find('output:') > -1): 
            f.write(f'```\n\n## {s_doc.strip()}\n```\n')
        elif (s_doc.find('description:') > -1): 
            f.write(f'```\n\n## {s_doc.strip()}\n```\n')
        else:
            f.write(f'{s_doc.strip()}\n')
    f.write('```')
    f.close()


# load data
mcds = pcdl.TimeStep(
    '../pcdl/data_timeseries_2d/output00000000.xml', 
    custom_type={'oncoprotein': str}, 
    verbose=False,
)

mcdsts = pcdl.TimeSeries(
    '../pcdl/data_timeseries_2d/', 
    custom_type={'oncoprotein': str}, 
    verbose=False,
)


# write makdown files
docstring_md(
    s_function = 'mcds.__init__', 
    ls_doc = pcdl.TimeStep.__init__.__doc__.split('\n'),
    s_header = "mcds = pcdl.TimeStep('path/to/outputnnnnnnnn.xml')"
)
docstring_md(
    s_function = 'mcds.set_verbose_false', 
    ls_doc = pcdl.TimeStep.set_verbose_false.__doc__.split('\n'), 
)
docstring_md(
    s_function = 'mcds.set_verbose_true', 
    ls_doc = pcdl.TimeStep.set_verbose_true.__doc__.split('\n'), 
)
docstring_md(
    s_function = 'mcds.get_multicellds_version', 
    ls_doc = pcdl.TimeStep.get_multicellds_version.__doc__.split('\n'), 
)
docstring_md(
    s_function = 'mcds.get_physicell_version', 
    ls_doc = pcdl.TimeStep.get_physicell_version.__doc__.split('\n'), 
)
docstring_md(
    s_function = 'mcds.get_timestamp', 
    ls_doc = pcdl.TimeStep.get_timestamp.__doc__.split('\n'), 
)
docstring_md(
    s_function = 'mcds.get_time', 
    ls_doc = pcdl.TimeStep.get_time.__doc__.split('\n'), 
)
docstring_md(
    s_function = 'mcds.get_runtime', 
    ls_doc = pcdl.TimeStep.get_runtime.__doc__.split('\n'), 
)
"""
docstring_md(
    s_function = 'mcds.', 
    ls_doc = pcdl.TimeStep..__doc__.split('\n'), 
)

docstring_md(
    s_function = 'mcds.', 
    ls_doc = pcdl.TimeStep..__doc__.split('\n'), 
)

docstring_md(
    s_function = 'mcds.', 
    ls_doc = pcdl.TimeStep..__doc__.split('\n'), 
)

docstring_md(
    s_function = 'mcds.', 
    ls_doc = pcdl.TimeStep..__doc__.split('\n'), 
)

docstring_md(
    s_function = 'mcds.', 
    ls_doc = pcdl.TimeStep..__doc__.split('\n'), 
)
"""

