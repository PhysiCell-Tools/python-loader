#os.system(f'python3 -c "import pcdl; help(pcdl.TimeStep.get_multicellds_version)" >> {s_opathfile}')


# library
import pcdl
import os


# function
def docstring_md(s_function, ls_doc, s_header=None, s_opath='./docstring/'):
    print(f'processing: {s_function} ...')
    os.makedirs(s_opath, exist_ok=True)
    s_opathfile = s_opath + f'{s_function}.md'
    f = open(s_opathfile, 'w')
    if (s_header is None):
        s_header = f'{s_function}()'
    f.write(f'# {s_header}\n\n')
    for s_doc in ls_doc:
        #if (len(s_doc.strip()) == 0):
        #    pass
        if (s_doc.find('input:') > -1):
            f.write(f'## {s_doc.strip()}\n```\n')
        elif (s_doc.find('output:') > -1):
            f.write(f'```\n\n## {s_doc.strip()}\n```\n')
        elif (s_doc.find('description:') > -1):
            f.write(f'```\n\n## {s_doc.strip()}\n```\n')
        else:
            f.write(f'{s_doc}\n')
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


# write pyMCDS initialize function makdown files
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


# write pyMCDS metadata function makdown files
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

# write pyMCDS mesh function markdown files
# range and axis
docstring_md(
    s_function = 'mcds.get_voxel_ijk_range',
    ls_doc = pcdl.TimeStep.get_voxel_ijk_range.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.get_mesh_mnp_range',
    ls_doc = pcdl.TimeStep.get_mesh_mnp_range.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.get_xyz_range',
    ls_doc = pcdl.TimeStep.get_xyz_range.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.get_voxel_ijk_axis',
    ls_doc = pcdl.TimeStep.get_voxel_ijk_axis.__doc__.split('\n'),
)
# mesh
docstring_md(
    s_function = 'mcds.get_mesh_mnp_axis',
    ls_doc = pcdl.TimeStep.get_mesh_mnp_axis.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.get_mesh',
    ls_doc = pcdl.TimeStep.get_mesh.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.get_mesh_2D',
    ls_doc = pcdl.TimeStep.get_mesh_2D.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.get_mesh_coordinate',
    ls_doc = pcdl.TimeStep.get_mesh_coordinate.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.get_mesh_spacing',
    ls_doc = pcdl.TimeStep.get_mesh_spacing.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.is_in_mesh',
    ls_doc = pcdl.TimeStep.is_in_mesh.__doc__.split('\n'),
)
# voxel
docstring_md(
    s_function = 'mcds.get_voxel_volume',
    ls_doc = pcdl.TimeStep.get_voxel_volume.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.get_voxel_spacing',
    ls_doc = pcdl.TimeStep.get_voxel_spacing.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.get_voxel_ijk',
    ls_doc = pcdl.TimeStep.get_voxel_ijk.__doc__.split('\n'),
)

# write pyMCDS microenv function markdown files
docstring_md(
    s_function = 'mcds.get_substrate_names',
    ls_doc = pcdl.TimeStep.get_substrate_names.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.get_substrate_dict',
    ls_doc = pcdl.TimeStep.get_substrate_dict.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.get_substrate_df',
    ls_doc = pcdl.TimeStep.get_substrate_df.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.get_concentration',
    ls_doc = pcdl.TimeStep.get_concentration.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.get_concentration_at',
    ls_doc = pcdl.TimeStep.get_concentration_at.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.get_conc_df',
    ls_doc = pcdl.TimeStep.get_conc_df.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.plot_contour',
    ls_doc = pcdl.TimeStep.plot_contour.__doc__.split('\n'),
)

# write pyMCDS cell agent function markdown files
docstring_md(
    s_function = 'mcds.get_cell_variables',
    ls_doc = pcdl.TimeStep.get_cell_variables.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.get_celltype_dict',
    ls_doc = pcdl.TimeStep.get_celltype_dict.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.get_cell_df',
    ls_doc = pcdl.TimeStep.get_cell_df.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.get_cell_df_at',
    ls_doc = pcdl.TimeStep.get_cell_df_at.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.plot_scatter',
    ls_doc = pcdl.TimeStep.plot_scatter.__doc__.split('\n'),
)

# write pyMCDS graph function markdown files
docstring_md(
    s_function = 'mcds.get_attached_graph_dict',
    ls_doc = pcdl.TimeStep.get_attached_graph_dict.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.get_neighbor_graph_dict',
    ls_doc = pcdl.TimeStep.get_neighbor_graph_dict.__doc__.split('\n'),
)
docstring_md(
    s_function = 'mcds.make_graph_gml',
    ls_doc = pcdl.TimeStep.make_graph_gml.__doc__.split('\n'),
)

# write pyMCDS unit function markdown files
docstring_md(
    s_function = 'mcds.get_unit_se',
    ls_doc = pcdl.TimeStep.get_unit_se.__doc__.split('\n'),
)


# write pyAnnData unit function markdown files
docstring_md(
    s_function = 'mcds.get_anndata',
    ls_doc = pcdl.TimeStep.get_anndata.__doc__.split('\n'),
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
