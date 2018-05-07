#script to parse a query file and show them in pymol
#by Matthew Baumgartner
#6-27-12

#usage:
#Please see the attached README file for detailed installation and usage instructions
#in pymol: 
#to install: `run /path/to/script/load_query.py`
#to run: `load_query <queryfile>`


#License:
#---------
#
#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.



#Change Log:
#-----------

#7-6-12 v0.1
# MB Prototype version, takes a .query file as an input and shows the spheres in the structure

#7-17-12 v0.1.1 MB
# Added support for getting the receptor and residues from the query file

#7-20-12 v0.1.2 MB
# Added support for using query files downloaded from ZincPharmer
# Draws arrows for the hydrogen bonds

#7-24-12 v0.1.3 MB
# Added better documentation
# Added command line args for whether to display the receptor and ligand and setting the transparency

#7-24-12 v0.1.4 MB
# If a feature is disabled in the query file, it will be disabled when it is loaded into pymol

#7-25-12 v0.1.5 MB
# Added loadall flag, which when set to False, does not load un-enabled points
# added string to boolean caster and better checking for whether the receptor is actually in the query file

#9-19-12 v1.0 MB
# Added support for running as a pymol plugin GUI
# Added option pane for controlling whether to change parameters (showing receptor and ligand, transparency etc.)

#2-21-13 v1.1 MB
# Added option for loading entire query as one object using 'as_one' option.
####TODO LIST
#Fix colors 7-6-12
#Plot arrows for HB acceptors and donors 7-6-12 - DONE 7-20-12 MB
#Test with downloaded query files from AQ/ZP 7-6-12 - DONE 7-20-12 MB
#load in the receptor from the query file if it exists as well 7-6-12 - DONE 7-20-12 MB
#write install script - checking os and making sure you can import simple json and adding to .pymolrc file 7-6-12
#add help the pymol can read - DONE 7-24-12
#check compatibility with anchorquery
#add option for loading the whole query as one object - DONE 2-21-13 MC
#load query as pseudo-atom so the mesh can be shown - prolly with arrows
#make work with the load command - if possible
#make menu option - DONE 8-25-12
#use auto_arg option
#load the anchor in nicely



from pymol import cmd
from pymol.cgo import *
import simplejson as json
import os
import urllib #for parsing the resides from the query file. Optional
import tempfile
import sys
import tkSimpleDialog
import tkMessageBox
import tkFileDialog
import Pmw
import Tkinter
import numpy as np
import math
import random


#a variable to determining how much output should be printed
#to the console
global DEBUG
DEBUG = 2



#convert strings to booleans
#after a bit of googleing, I can't find a built in python operator to do this, replace this if there is
def str_to_bool(string):
    """
    Convert an input string to a boolean
    
    Arguments:
    string -- string to be cast to boolean ('True' or 'False')
    
    Returns:
    boolean -- boolean representation of string
    """
    #if it's already a boolean, return it
    if isinstance(string, bool):
        return string
    #if its a string, but not True or False, throw an error
    elif string not in ['True','False']:
        sys.stderr.write('Error: input parameter: ' + string + " must be 'True' or 'False'\n")
        sys.exit()
    else:
        if string == 'True':
            return True
        elif string == 'False':
            return False
        

#Draw a arrow for the pharmacophore
#partially adapted from: cones_native.py -> (c)2006 Tsjerk A. Wassenaar, PhD, University of Utrecht
#Available here: rna.ucsc.edu/rnabasepairs/cones_native.py as of 7-17-2012 MB
def draw_pharma_arrow(centx, centy, centz, unitx, unity, unitz, radius, startcolor = [1,1,1], endcolor = [1,1,1], cyllen = 1.0, cylrad = 0.4, conelen = 0.5, conerad = 0.5):
    '''
    Create a pymol cgo arrow object made from a cylinder and a cone.
    
    Arguments:
    centx, centy, centz -- coordinates of center of pharmacophore
    unitx, unity, unitz -- unit vector - length = 1
    radius -- the radius of the pharmacophore sphere
    
    Keyword Arguments:
    startcolor -- starting color for gradient coloring the cone (rgb)
    endcolor -- starting color for gradient coloring the cone (rgb)
    cyllen -- length of the cylinder portion of the arrow (default = 1.0)
    cylrad -- radius of the cylinder portion of the arrow (default = 0.5)
    conelen -- length of the cone portion of the arrow (default = 0.5)
    conerad -- radius of the base of the cone (default = 0.5)
    
    Returns:
    arrowlist -- list of floats that specify the arrow to be drawn
    '''
    
    
    #get the point on the vector that is on the radius
    #make the radius a bit smaller so the graphics line up better 
    basex = centx + unitx * (radius * 0.90)
    basey = centy + unity * (radius * 0.90)
    basez = centz + unitz * (radius * 0.90)
    
    #get the point that is the end of the cylinder and is the base of the cone
    midx = basex + unitx * cyllen
    midy = basey + unity * cyllen
    midz = basez + unitz * cyllen
    
    #get the point that is the end of the cone
    endx = midx + unitx * conelen
    endy = midy + unity * conelen
    endz = midz + unitz * conelen
    
    arrowlist = [CYLINDER, basex, basey, basez, midx, midy, midz, cylrad] + startcolor + endcolor
    #                                               cone radii:  base     tip
    #                  type  start x, y,   z    end x,  y,   z,    ^       ^       start and end colors   ? state?
    arrowlist.extend([ CONE , midx, midy, midz, endx, endy, endz, conerad, 0] + startcolor + endcolor + [ 1,1 ])

    return arrowlist

def compute_normal(x1, y1, z1, x2, y2, z2, x3, y3, z3):
    '''
    # A helper function for computing the normal to a triangular facet
    by Gareth Stockwell, 2004
    
    Found here: http://cavanagh-lab.bch.ncsu.edu/new/intranet/pymol/box.py
    '''
    nx = (y2-y1)*(z3-z2) - (z2-z1)*(y3-y2)
    ny = (z2-z1)*(x3-x2) - (x2-x1)*(z3-z2)
    nz = (x2-x1)*(y3-y2) - (y2-y1)*(x3-x2)

    return (nx,ny,nz)

def draw_triangle(x1,y1,z1, x2,y2,z2, x3,y3,z3, c1,c2,c3):
    '''
    Generate list for drawing a cgo triangle.
    
    Arguments:
    x1,y1,z1 -- coordinates of the first vertex
    x2,y2,z2 -- coordinates of the second vertex
    x3,y3,z3 -- coordinates of the third vertex
    c1,c2,c3 -- RGB color values
    
    Returns:
    trianglelist -- list of cgo floats
    '''
    
    #compute the normals
    n1,n2,n3 = compute_normal(x1, y1, z1, x2, y2, z2, x3, y3, z3)
    
#    trianglelist = [BEGIN, TRIANGLES, NORMAL, n1,n2,n3, VERTEX, x1,y1,z1, VERTEX, x2,y2,z2, VERTEX, x3,y3,z3, END]
    trianglelist = [BEGIN, TRIANGLES, COLOR, c1,c2,c3, VERTEX, x1,y1,z1, VERTEX, x2,y2,z2, VERTEX, x3,y3,z3, END]
    
    return trianglelist

def generate_sphere_points(n):
    """
    Returns list of 3d coordinates of points on a sphere using the
    Golden Section Spiral algorithm.
    
    Arguments:
    n -- number of points to generate
    
    Returns:
    points -- np array of the points (n x 3)
    """
    points = []
    inc = math.pi * (3 - math.sqrt(5))
    offset = 2 / float(n)
    for k in range(int(n)):
        y = k * offset - 1 + (offset / 2)
        r = math.sqrt(1 - y*y)
        phi = k * inc
        points.append([math.cos(phi)*r, y, math.sin(phi)*r])
    return np.array(points) #MB

def get_neighbors(points, point, num_neighbors):
    '''
    Naive method for finding the nearest neighbors of a
    point on a sphere of points
    
    Arguments:
    points -- list of points (numpy arrays)
    point -- point of interest
    num_neighbors -- number of neighbors to return
    
    Returns:
    neighbors -- list of neighboring points
    '''

    distances = [ [pt, np.linalg.norm(point - pt)] for pt in points]
    
    distances = sorted(distances, key = lambda foo: foo[1])

    #don't return the point itself - i.e. the first one in the list
    neighbors = np.array(distances)[1:(num_neighbors + 1), 0]

    return neighbors
    
def make_mesh(x,y,z, rad, color = [0,0,0], n = 100, num_neighbors = 6, rand_rot = True):
    '''
    Draw a wire mesh sphere of radius rad centered at x,y,z
    
    Arguments:
    x,y,z -- coordinates of center of sphere
    rad -- radius of the sphere
    n -- number of sphere points -> density of the mesh
    neighbors -- number of nearest neighbors to draw lines with
    rand_rot -- boolean for whether to rotate the spheres a little.
                This allows for viewing of overlapping pharmacophores

    Returns:
    mesh -- list of cgo floats
    '''
    #get the points centered at the origen
    points = generate_sphere_points(n)
    
    #expand based on the radius
    points = points * rad
    
    if rand_rot:
        #generate a random amount to rotate (0 to 2*pi) around the z axis
        rot = random.random() * 2*math.pi
        
        #set up the transformation matrix
        trans = np.array([[math.cos(rot), -math.sin(rot), 0, 0], \
                          [math.sin(rot), math.cos(rot), 0, 0], \
                          [0, 0, 1, 0], \
                          [0, 0, 0, 1] ])
        
        #append a 1 to each coordinate vector for doing the transformation
        tmppoints = [ np.array([ pt[0], pt[1], pt[2], 1]) for pt in points ]
        
        
        #apply the rotation
        tmppoints = np.array([ trans.dot(pt) for pt in tmppoints ])
        
        #Strip off the 1 that was added
        points = tmppoints[:,:3] 
        
    #center them at the specified location
    points = points + np.array([x,y,z])
        
    #start building the cgo list
    mesh = [LINES]
    mesh.extend([COLOR, color[0], color[1], color[2]])
    
    for point in points:
        neighbors = get_neighbors(points, point, num_neighbors)
        for nbr in neighbors:
            mesh.extend([VERTEX, point[0], point[1], point[2]])
            mesh.extend([VERTEX, nbr[0], nbr[1], nbr[2]])
        
    return mesh

def load_query(query_file = '', receptor = True, ligand = True, trans = 0.7, loadall = True, mesh = True, as_one = False):
    '''
DESCRIPTION
    
    Read a .query json-encoded session file from ZincPharmer (http://zincpharmer.csb.pitt.edu/), and display the pharmacophore spheres.

USAGE

    load_query <query_file> [, receptor = True [, ligand = True [, loadall = True[, trans = 0.7 [, mesh = True]]]]

ARGUMENTS
    query_file -- The zincpharmer session file (*.query)
    receptor -- By default, if there is a receptor in the query file, it will be loaded. Set to False to disable
    ligand -- By default, if there is a ligand in the query file, it will be loaded. Set to False to disable
    trans -- Transparency of the pharmacophore spheres (0.0 - 1.0)
    loadall -- If False, only load points that are enabled, else load all (Default: True)
    mesh -- If True, make a wire mesh to represent the pharmacophore
    as_one -- If True, load the entire set of pharmacophore features as one object. Useful for comparing two queries. Overrides loadall, receptor and ligand.
    
NOTES
    Online example query file: http://zincpharmer.csb.pitt.edu/examples/kinase.query
    See the included Readme file for more info
    
EXAMPLES
    load_query pharmer.query
    load_query /path/to/file/kinase.query, receptor = False
    
PYTHON API
    load_query(query_file, receptor = 'True', ligand = 'True', trans = 0.7, loadall = True, as_one = False):
    '''
    
    RECEPTOR = str_to_bool(receptor) #boolean for whether to load the receptor
    LIGAND = str_to_bool(ligand) #boolean for whether to load the ligand
    TRANS = float(trans) #Transparency of the spheres
    LOADALL = str_to_bool(loadall) #whether to load all the points or just the enabled ones
    USE_MESH = str_to_bool(mesh)
    AS_ONE = str_to_bool(as_one)
    
    #if you are loading all the points as one object, turn off loading all points as well as the receptor and ligand
    if AS_ONE:
        LOADALL = False
        RECEPTOR = False
        LIGAND = False
    
    #parameters for drawing the arrows - these draw pretty good arrows, but if you need to change them do it here
    cyllen = 0.70
    cylrad = 0.1
    conelen = 0.5
    conerad = 0.25
    
    #RGB colors for coloring the pharmacophores
    #'Other' is default - if unknown type
    colorlist = { 'Hydrophobic':        [0.100, 1.000, 0.000],
                  'HydrogenAcceptor':   [1.000, 0.84, 0.000],\
                  'HydrogenDonor':      [1.000, 1.000, 1.000],\
                  'Aromatic':           [0.627, 0.1254, 0.941],\
                  'NegativeIon':        [1.00, 0.00, 0.00],\
                  'PositiveIon':        [0.00, 0.00, 1.00],\
                  'Other':              [0.74, 0.74, 0.74],\
                  'PhenylalanineAnalog':[1.0, 1.0, 0.0],\
                  'LeuValAnalog':       [1.0, 1.0, 0.0] \
                  }


    #Edit some pymol settings
    #set the gui width to be bigger so the object names don't get truncated
    cmd.set('internal_gui_width', '400')
    
    #load in the query file
    qfile = open(query_file).read()
    #parse it with json
    pts = json.loads(qfile)
    
    #print 'RECEPTOR:',RECEPTOR
    
    if RECEPTOR:
        #if there is a receptor record, pull it out
        if pts.has_key('receptor'):
            
            #receptor pdb as one string
            rec = pts.get('receptor')
    
            #make sure the receptor record is not nulled
            if rec not in ['', None]:            
                #open a new temp file and write the receptor to it 
                tfout = tempfile.NamedTemporaryFile(suffix = '.pdb', delete=False)
                tmpname = tfout.name 
                tfout.write(rec)
                tfout.close()
    #            print 'Receptor file saved:', tmpname
        
                #Try and get the recname from the query file else default to the base name of the input file       
                if pts.has_key('recname'):
                    recfname = pts.get('recname')
                else:
                    recfname = os.path.basename(file)
                
                #load the receptor file in
                cmd.load(tmpname, recfname)
                
                #show the surface on the receptor
                cmd.show('surface', recfname)
                
                #check to see if the temp file exists and if so delete it
                if DEBUG > 2:
                    print 'tfout:', tfout.name
                if os.path.exists(tfout.name):
                    os.remove(tfout.name)

    if LIGAND:
        #load in the residues if they exist
        if pts.has_key('sdf'):
            if pts.get('sdf') not in  ['',None]:
                sdf = pts.get('sdf')
                
                #decode the url
                sdf = urllib.unquote(sdf) 
                #write it out to a temp file
                tfout = tempfile.NamedTemporaryFile(suffix = '.pdb', delete=False)
                tmpname = tfout.name 
                tfout.write(sdf)
                tfout.close()
    
                #load it into pymol
                cmd.load(tmpname, 'residues')
                #show as sticks
                cmd.show('sticks', 'residues')
                #hide the hydrogens on the residues
                cmd.hide('everything', 'symbol H and /residues')
            
                cmd.orient('residues')
                
                #check to see if the temp file exists and if so delete it
                if os.path.exists(tfout.name):
                    os.remove(tfout.name)
    
    #list to hold all the query points if using AS_ONE
    whole_query = []
    
    #if there is an anchor in the query file, try to display it
    if pts.has_key('anchors'):
        for anchor in pts.get('anchors'):
            name  = anchor.get('name')
            #not sure what to do with this yet, but grab it for now
            if anchor.has_key('wiggle'):
                wiggle = anchor.get('wiggle')
            
            enab = True
            if anchor.has_key('enabled'):
                enab = anchor.get('enabled') 
            
            x1 = anchor.get('x1')
            y1 = anchor.get('y1')
            z1 = anchor.get('z1')
            
            x2 = anchor.get('x2')
            y2 = anchor.get('y2')
            z2 = anchor.get('z2')
            
            x3 = anchor.get('x3')
            y3 = anchor.get('y3')
            z3 = anchor.get('z3')
            
            c1,c2,c3 = colorlist[name]
            
            tri = draw_triangle(x1,y1,z1, x2,y2,z2, x3,y3,z3, c1,c2,c3)
    
#            print 'BEGIN, TRIANGLE, COLOR, c1,c2,c3, NORMAL, n1,n2,n3, VERTEX, x1,y1,z1, VERTEX, x2,y2,z2, VERTEX, x3,y3,z3, END'
#            print tri
#            sys.stderr.write(','.join(tri) + '\n')
    
    
            if AS_ONE:
                if enab:
                    whole_query.append(tri)
            else:
                cmd.load_cgo(tri, name)
                if not enab:
                    cmd.disable(name)
                
                
    #for each point, pull out the coordinates and radius
    if pts.has_key('points'):
        for p in pts.get('points'):
            name = p.get('name')
            
            x,y,z = p.get('x'), p.get('y'), p.get('z')
            rad = p.get('radius')
            
            #get the enabled state
            enab = True
            if p.has_key('enabled'):
                enab = p.get('enabled')
            
            if DEBUG > 0:
                print name, x,y,z, 'enabled =', enab, 'rad =', rad
            
            if LOADALL or enab: 
                #make sure that the pharmacophore type is a known one, else color it Grey
                if name in colorlist:
                    colors = colorlist[name]
                else:
                    colors = colorlist['Other']
    
                #name that will be displayed in pymol's internal gui menu            
                dispname = name + '_' + str(round(x,3)) + '_' + str(round(y,3)) + '_' + str(round(z,3))
                
                if USE_MESH:
                    spherelist = [BEGIN]
                    spherelist += make_mesh(x,y,z, rad, color = colors)
                else:
                    #Build list for the pharmacophore sphere
                    spherelist = [BEGIN, ALPHA, TRANS, COLOR] + colors + [SPHERE, x, y, z, rad]
                
                if DEBUG > 4:
                    print '\n\n\n\n\n\n\nspherelist:', spherelist
                
                #clear the arrow list
                arrowlist = []
                
                #plot arrows if there are vectors
                #The arrow will be plotted from the surface of the radius, out to the end of the line of the arrow for a 
                #specified distance, then a cone at the end of that. We are given the center and a vector of length 1 in the 
                #direction. it will look like this eventually
                #Given vector with radius as ')':
                #  .----->  )    
                #want to draw an arrow like
                #  .        )----->
                if p.has_key('svector'):
                    svec = p.get('svector')
                    if svec != None:
                        #input vector has length of 1
                        svx = svec.get('x')
                        svy = svec.get('y')
                        svz = svec.get('z')
                        
                        arrowlist = draw_pharma_arrow(x,y,z, svx,svy,svz, rad, startcolor=colors, endcolor=colors, cyllen=cyllen, cylrad=cylrad, conelen=conelen, conerad=conerad )
                
                
                if AS_ONE:
                    if enab:
                        whole_query.extend(spherelist + arrowlist)
                else:
                    cmd.load_cgo(spherelist + arrowlist, dispname, 1)
                    if not enab:
                        cmd.disable(dispname)
    if AS_ONE:
        if DEBUG > 4:
            print 'whole_query'
            print whole_query
        cmd.load_cgo(whole_query, os.path.basename(query_file),1)
           
    #refresh the plot so you see everything
    cmd.refresh()

cmd.extend("load_query", load_query)

default_settings = {
            'receptor' : True,
            'ligand' : True,
            'transparency' : 0.7,
            'asone' : False
            }

class LoadQuery:
    ''' Main plugin class '''
    def __init__(self, app):
        parent = app.root
        self.parent = parent
        
        self.query_file = self.load_query_dialog(app)
        
        if not self.query_file:
            return
        
        self.show_receptor = default_settings['receptor']
        self.show_ligand = default_settings['ligand']
        self.asone = default_settings['asone']
        
        #need to use a tkinter.doublevar object to store the transparency so 
        #that it will be displayed in the GUI 
        self.transparency = Tkinter.DoubleVar()
        self.transparency.set(default_settings['transparency'])
        
        
        self.display_options = {
                'receptor' : default_settings['receptor'],
                'ligand' : default_settings['ligand'],
                'transparency' : default_settings['transparency'],
                'asone': default_settings['asone']
                }
        
        self.dialog = Pmw.Dialog(parent, 
                                 buttons = ('Ok', 'Exit'), 
                                 title = 'Pharmacophore Query Loader',
                                 command = self.button_pressed )
    
        self.dialog.withdraw()
    
        Pmw.setbusycursorattributes(self.dialog.component('hull'))
    
        self.dialog.geometry('650x480')
        self.dialog.bind('<Return>',self.button_pressed)
    
        title_label = Tkinter.Label(self.dialog.interior(),
                                             text = 'PyMOL Pharmacophore Query Viewer\nMatthew Baumgartner\nmpb21@pitt.edu',
                                             background = 'navy',
                                             foreground = 'white',
                                             )
    
        title_label.pack(expand = 0, fill = 'both', padx = 4, pady = 4)
    
        #main notebook
        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill = 'both', expand = 1, padx = 3, pady = 3)
    
        #build option page
        self.option_page = self.notebook.add('Options')
        
        self.option_page_display = Pmw.Group(self.option_page, tag_text='Options')
        self.option_page_display.pack(fill = 'x', expand = 1, padx = 10, pady = 0)
        
        
        self.load_receptor_radio = Pmw.RadioSelect(self.option_page_display.interior(),
                                         selectmode='multiple',
                                         buttontype='checkbutton',
                                         labelpos='w',
                                         label_text='Display mode',
                                         orient='vertical',
                                         frame_relief='ridge',
                                         command=self.load_receptor_mode_changed)
        self.load_receptor_radio.pack(side = 'top', padx = 10, anchor='w')
        
        
        for entry, label in (('asone', 'Load as one object'), ('receptor','Load Receptor'), ('ligand', 'Load Ligand')):
            self.load_receptor_radio.add(entry, text = label)
        for entry in ('asone', 'receptor','ligand'):
            if self.display_options[entry]:
                self.load_receptor_radio.invoke(entry)
                
        if self.display_options['asone']:
            self.load_receptor_radio.deselect('receptor')
            self.load_receptor_radio.deselect('ligand')
                
        #Add a input for the transparency
        
        self.transparency_frame = Tkinter.Frame(self.option_page_display.interior())
        self.transparency_label = Tkinter.Label(self.transparency_frame, text = 'Transparency')
        self.transparency_location = Tkinter.Entry(self.transparency_frame, 
                                                   textvariable = self.transparency, 
                                                   bg = 'black', 
                                                   fg = 'white', 
                                                   width = 10)
        self.transparency_scrollbar = Tkinter.Scrollbar(self.transparency_frame, orient='horizontal', command = self.transparency_changed)
        
        self.transparency_label.pack(side = 'left')
        self.transparency_location.pack(side = 'left', padx = 5)
        self.transparency_scrollbar.pack(side = 'left')
        self.transparency_frame.pack(side = 'left', padx = 4, pady = 1)
        
        
        
        
#        Pmw.alignlabels([self.load_receptor_radio, self.transparency_frame])
#        Pmw.aligngrouptags(self.option_page_display)
        
        #done
#        parent.mainloop()
        self.notebook.setnaturalsize()
        self.dialog.show()
        
        if DEBUG > 3:
            print 'display_options:'
            for key,val in self.display_options.iteritems():
                print key, val
        
        
    
    def load_query_dialog(self, app):
        file = tkFileDialog.askopenfilename(title = 'Open query file...',
                                            defaultextension = '.query',
                                            filetypes = [('Query files', '.query'), ('All Files', '*.*')]
                                            )
        if file:
            #open the dialog box
            #self.load_options_dialog(app)
            
            return file
        else:
            return None
            
    def button_pressed(self, result):
        if DEBUG > 2:
            print 'button_pressed result', result
        if hasattr(result,'keycode'):
            if result.keycode == 36:
                print 'keycode:', result.keycode
        elif result == 'Ok':
            if DEBUG > 2:
                print 'loading query'
                print self.query_file, 'receptor =', self.display_options['receptor'], 'ligand =', self.display_options['ligand']
            load_query(self.query_file, 
                        receptor = self.display_options['receptor'],
                        ligand = self.display_options['ligand'],
                        trans = self.transparency.get(),
                        as_one = self.display_options['asone'])
            self.dialog.withdraw()
                        
        elif result == 'Exit' or result == None:
            self.dialog.withdraw()

                
    def load_receptor_mode_changed(self, button_name, pressed):
        if pressed:
            self.display_options[button_name] = True
            if button_name == 'asone':
                #as_one overrides receptor and ligand
                self.display_options['receptor'] = False
                self.display_options['ligand'] = False
                
                #deselect the receptor and ligand buttons
                if 'receptor' in self.load_receptor_radio.getvalue():
                    self.load_receptor_radio.invoke('receptor')
                if 'ligand' in self.load_receptor_radio.getvalue():
                    self.load_receptor_radio.invoke('ligand')
                    
            action = 'Enabled'
        else:
            self.display_options[button_name] = False
            action = 'Disabled'
        txt = action+' ligand display mode <'+button_name+'>'
#        status_line.configure(text=txt)
        if DEBUG > 2:
            print txt
        
    def transparency_changed(self, x):
        incr = 0.05
#        old_trans = self.display_options['transparency']
        old_trans = self.transparency.get()
        
        new_trans = old_trans + float(x)*incr
        #round it off to 2 digits
        new_trans = round(new_trans, 2)
        if new_trans >= 1.0:
            new_trans = 1.0
        elif new_trans <= 0.0:
            new_trans = 0.0
        if DEBUG > 2:
            print 'changing transparency from:', old_trans, 'to:', new_trans
        
#        self.display_options['transparency'] = new_trans
        self.transparency.set(new_trans)
        
        
def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command', 'Pharmacophore Query Loader', label = 'Pharmacophore Query Loader', command = lambda s=self : LoadQuery(s))
    

#From http://pymol.sourceforge.net/newman/user/S0500cgo.html
#
#A CGO is simply a Python list of floating point numbers, which are compiled by PyMOL into a CGO object and associated with a given state.
#
#Lowercase names below are should be replaced with floating-point numbers. Generally, the TRIANGLE primitive should only be used only as a last restore since it is much less effective to render than using a series of vertices with a BEGIN/END group.
#
#BEGIN, { POINTS | LINES | LINE_LOOP | LINE_STRIP | TRIANGLES | TRIANGLE_STRIP | TRIANGLE_FAN },
#
#VERTEX, x,  y,  z,
#
#COLOR,  red, green, blue, 
#
#NORMAL, normal-x, normal-y,  normal-z, 
#
#END,
#
#LINEWIDTH, line-width, 
#
#WIDTHSCALE, width-scale,   # for ray-tracing
#
#SPHERE, x, y, z,  radius    # uses the current color
#
#CYLINDER, x1, y1, z1, x2, y2, z2, radius,
#          red1, green1, blue1, red2, green2, blue2,
#
#TRIANGLE,  x1, y1, z1, 
#           x2, y2, z2,
#           x3, y3, z3,
#           normal-x1, normal-y1, normal-z1,
#           normal-x2, normal-y2, normal-z2,
#           normal-x3, normal-y3, normal-z3,
#           red1, green1, blue1,          
#           red2, green2, blue2,          
#           red3, green3, blue3,          
