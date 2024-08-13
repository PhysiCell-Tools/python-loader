#### MCDS Time Steps, Pandas, and Plotting                                      
                                                                                
Since microenvironment data and cell data can be retrieved as pandas datafarme, basic plotting (line plot, bar plot, histogram, boxplot, kernel density estimation plot, area plot, pie plot, scatter plot, hexbin plot) can easily be generated with the **pandas plot function**.
As mentioned above, for microenvironment data, the TimeStep class has a mcds.get\_contour function because pandas has no contour and contourf plots implementation.
All these plots are **mathplotlib** plots, hence fine tuning can always be done using the matplotlib library.
                                                                                
+ https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.plot.html
+ https://matplotlib.org/                                                       
                                                                                
```python                                                                       
# load library                                                                  
import matplotlib.pyplot as plt                                                 
                                                                                
# set constantes                                                                
z_slice = 0                                                                     
                                                                                
# generate plot                                                                 
fig, ax = plt.subplots(figsize=(7,4))                                           
                                                                                
# plot microenvironment                                                         
mcds.plot_contour('oxygen', z_slice=z_slice, ax=ax)                             
                                                                                
# plot cells                                                                    
df = mcds.get_cell_df()                                                         
df = df.loc[(df.mesh_center_p == z_slice),:]                                    
df.plot(                                                                        
    kind = 'scatter',                                                           
    x = 'position_x',                                                           
    y = 'position_y',                                                           
    s = 10,                                                                     
    c = 'oncoprotein',                                                          
    cmap = 'magma',                                                             
    vmin = 0,                                                                   
    vmax = 2,                                                                   
    grid = True,                                                                
    title = 'cells and microenvironment',                                       
    ax = ax,                                                                    
)                                                                               
ax.axis('equal')                                                                
plt.tight_layout()                                                              
                                                                                
# save plot to file                                                             
# bue 20230101: note, in this test dataset, cells were seeded outside the actual domain.
fig.savefig('pymcds_2d_cell_and_microenvironment.png', facecolor='white')       
plt.close()                                                                        
```                                                                             

