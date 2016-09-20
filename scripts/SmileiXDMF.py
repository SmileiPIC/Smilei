import h5py
import numpy as np
import os

from jinja2 import Template

def fieldsXDMF():

    field2D = Template("""
    <Xdmf Version="2.0">
        <Domain Name="{{field}}">
            <Grid GridType="Collection" CollectionType="Temporal">
                {% for step,time in steps_times -%}
                <Grid Name="Mesh" GridType="Uniform">
                    <Time Value="{{time}}" />
                    <Topology TopologyType="2dCorectMesh" Dimensions="{{size[0]}} {{size[1]}}"/>
                    <Geometry GeometryType="ORIGIN_DXDY">
                        <DataItem Format="XML" NumberType="float" Dimensions="2">{{orig[0]}} {{orig[1]}}</DataItem>
                        <DataItem Format="XML" NumberType="float" Dimensions="2">{{cell_length[0]}} {{cell_length[1]}}</DataItem>
                    </Geometry>
                    <Attribute Name="{{field}}" Center="Node" AttributeType="Scalar">
                        <DataItem NumberType="Float" Precision="8" Dimensions="{{size[0]}} {{size[1]}}" Format="HDF">../{{filename}}:/{{step}}/{{field}}</DataItem>
                    </Attribute>
                </Grid>
                {% endfor -%}
            </Grid>
        </Domain>
    </Xdmf>""")
    field3D = Template("""
    <Xdmf Version="2.0">
        <Domain Name="{{field}}">
            <Grid GridType="Collection" CollectionType="Temporal">
                {% for step,time in steps_times -%}
                <Grid Name="Mesh" GridType="Uniform">
                    <Time Value="{{time}}" />
                    <Topology TopologyType="3dCorectMesh" Dimensions="{{size[0]}} {{size[1]}} {{size[2]}}"/>
                    <Geometry GeometryType="ORIGIN_DXDYDZ">
                        <DataItem Format="XML" NumberType="float" Dimensions="3">{{orig[0]}} {{orig[1]}} {{orig[2]}}</DataItem>
                        <DataItem Format="XML" NumberType="float" Dimensions="3">{{cell_length[0]}} {{cell_length[1]}} {{cell_length[2]}}</DataItem>
                    </Geometry>
                    <Attribute Name="{{field}}" Center="Node" AttributeType="Scalar">
                        <DataItem NumberType="Float" Precision="8" Dimensions="{{size[0]}} {{size[1]}} {{size[2]}}" Format="HDF">../{{filename}}:/{{step}}/{{field}}</DataItem>
                    </Attribute>
                </Grid>
                {% endfor -%}
            </Grid>
        </Domain>
    </Xdmf>""")

    for fname in ['Fields.h5','Fields_avg.h5']:
        print "File ", fname
        if os.path.exists(fname):
    
            fh5 = h5py.File(fname, 'r')

            fields = fh5.itervalues().next().keys() # list of fields

            # array extract from fh5.keys() except last (tmp)
            keys = fh5.keys()[0:len(fh5.keys())-1]
            timesteps = np.array(keys)

            my_steps_times=zip(timesteps,timesteps.astype(np.float)*fh5.attrs['res_time'])

            celllength=fh5.attrs['cell_length']

            if len(celllength) == 2:
                my_tmpl=field2D
            elif len(celllength) == 3:
                my_tmpl=field3D
            else:
                raise "Error dimension not supported"

            for myfield in fields:
                print "\tField", myfield

                sizes=fh5[fh5.keys()[0]+'/'+myfield].shape
                if myfield.startswith('Ex') :
                    origin=[-0.5,0,0]
                elif myfield.startswith('Ey') :
                    origin=[0,-0.5,0]
                elif myfield.startswith('Ez') :
                    origin=[0,0,-0.5]
                elif myfield.startswith('Bx') :
                    origin=[0,-0.5,-0.5]
                elif myfield.startswith('By') :
                    origin=[-0.5,0,-0.5]
                elif myfield.startswith('Bz') :
                    origin=[-0.5,-0.5,0]
                elif myfield.startswith('Jx') :
                    origin=[-0.5,0,0]
                elif myfield.startswith('Jy') :
                    origin=[0,-0.5,0]
                elif myfield.startswith('Jz') :
                    origin=[0,0,-0.5]
                elif myfield.startswith('Rho') :
                    origin=[0,0,0]
                else :
                    raise "This should not happend: implement",myfield
            
                origin=origin[:len(celllength)]
                with open(out_dir+'/'+os.path.splitext(os.path.basename(fname))[0]+'_'+myfield+'.xdmf','w') as fout:

                    render=my_tmpl.render(
                        filename=fh5.file.filename,
                        steps_times=my_steps_times,
                        orig=origin*celllength,
                        field=myfield,
                        size=sizes,
                        cell_length=celllength,
                        sim_length=fh5.attrs['sim_length'])
    
                    fout.write(render)
                    fout.close()

def probesXDMF():
    fh5 = h5py.File('Probes.h5', 'r')
    
    
    probe2D = Template("""
    <Xdmf Version="2.0">
        <Domain>
            <Grid GridType="Collection" CollectionType="Temporal">
                {% for step,time in steps_times -%}
                <Grid Name="Mesh" GridType="Uniform">
                    <Time Value="{{time}}"/>
                    <Topology TopologyType="2dSMesh" Dimensions="{{size[1]}} {{size[0]}}"/>
                    <Geometry GeometryType="XY">
                        <DataItem Format="XML" NumberType="float" Dimensions="{{size[0]*size[1]}} 2">
                        {% for i in range(prob[1].get("positions").shape[1]) -%}{{prob[1].get("positions").value[0][i]}} {{prob[1].get("positions").value[1][i]}} {% endfor %}
                        </DataItem>
                    </Geometry>
                    {% for i in range(10) -%}
                    <Attribute Name="{{conv[i]}}" Center="Node" AttributeType="Scalar">
                        <DataItem ItemType="HyperSlab" Dimensions="{{size[0]*size[1]}}" Type="HyperSlab">
                            <DataItem Dimensions="3 2" Format="XML">
                                0 {{i}}
                                1  1
                                {{size[0]*size[1]}} 1
                            </DataItem>
                            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="{{size[0]*size[1]}} 10">
                                ./Probes.h5:/{{prob[0]}}/{{step}}
                            </DataItem>
                        </DataItem>
                    </Attribute>
                    {% endfor -%}
                </Grid>
                {% endfor -%}
            </Grid>
        </Domain>
    </Xdmf>""")

# <!--                         <DataItem NumberType="Float" Precision="8" Dimensions="10 {{size[1]}} {{size[0]}}" Format="HDF">../Probes.h5:/{{probe_name}}/{{step}}</DataItem> -->
    probe3D = Template("""""")
#     """
#     <Xdmf Version="2.0">
#         <Domain Name="{{field}}">
#             <Grid GridType="Collection" CollectionType="Temporal">
#                 {% for step,time in steps_times -%}
#                 <Grid Name="Mesh" GridType="Uniform">
#                     <Time Value="{{time}}" />
#                     <Topology TopologyType="3dCorectMesh" Dimensions="{{size[0]}} {{size[1]}} {{size[2]}}"/>
#                     <Geometry GeometryType="ORIGIN_DXDYDZ">
#                         <DataItem Format="XML" NumberType="float" Dimensions="3">{{orig[0]}} {{orig[1]}} {{orig[2]}}</DataItem>
#                         <DataItem Format="XML" NumberType="float" Dimensions="3">{{cell_length[0]}} {{cell_length[1]}} {{cell_length[2]}}</DataItem>
#                     </Geometry>
#                     <Attribute Name="{{field}}" Center="Node" AttributeType="Scalar">
#                         <DataItem NumberType="Float" Precision="8" Dimensions="{{size[0]}} {{size[1]}} {{size[2]}}" Format="HDF">../{{filename}}:/{{step}}/{{field}}</DataItem>
#                     </Attribute>
#                 </Grid>
#                 {% endfor -%}
#             </Grid>
#         </Domain>
#     </Xdmf>""")

    for probe in fh5.items():
        fname=probe[0]

        my_shape=probe[1].get("number").value
        
        print "Probe", fname, my_shape

        if len(my_shape) == 2:
            my_tmpl=probe2D
        elif len(my_shape) == 3:
            my_tmpl=probe3D
        else:
            raise "Error dimension not supported"

        times = []
        timesstr = []
        for key in probe[1].iterkeys():
            try   : 
                times.append( int(key) )
                timesstr.append(key)
            except: pass

        my_steps_times=zip(timesstr,times)
        
        with open(out_dir+'/Probe_'+probe[0]+'.xdmf','w') as fout:
            render=my_tmpl.render(
                prob=probe,
                steps_times=my_steps_times,
                size=my_shape,
                conv = ["Ex","Ey","Ez","Bx","By","Bz","Jx","Jy","Jz","Rho"])

            fout.write(render)
            fout.close()
                    

    print "Probes done"

##############################################################
out_dir='smileiXDMF'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)





#probesXDMF()

fieldsXDMF()


