#!/usr/bin/env python2

"""Python script to plot tally data generated by OpenMC."""

import os
import sys

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np

from statepoint import *

if sys.version_info[0] < 3:
    import Tkinter as tk
else:
    import tkinter as tk
import tkFileDialog
import tkFont
import tkMessageBox
import ttk


class MeshPlotter(tk.Frame):
    def __init__(self, parent, filename):
        tk.Frame.__init__(self, parent)

        self.labels = {'cell': 'Cell:', 'cellborn': 'Cell born:',
                       'surface': 'Surface:', 'material': 'Material:',
                       'universe': 'Universe:', 'energyin': 'Energy in:',
                       'energyout': 'Energy out:'}

        self.filterBoxes = {}

        # Read data from source or leakage fraction file
        self.get_file_data(filename)

        # Set up top-level window
        top = self.winfo_toplevel()
        top.title('Mesh Tally Plotter: ' + filename)
        top.rowconfigure(0, weight=1)
        top.columnconfigure(0, weight=1)
        self.grid(sticky=tk.W+tk.N)

        # Create widgets and draw to screen
        self.create_widgets()
        self.update()

    def create_widgets(self):
        figureFrame = tk.Frame(self)
        figureFrame.grid(row=0, column=0)

        # Create the Figure and Canvas
        self.dpi = 100
        self.fig = Figure((5.0, 5.0), dpi=self.dpi)
        self.canvas = FigureCanvasTkAgg(self.fig, master=figureFrame)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Create the navigation toolbar, tied to the canvas
        self.mpl_toolbar = NavigationToolbar2TkAgg(self.canvas, figureFrame)
        self.mpl_toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Create frame for comboboxes
        self.selectFrame = tk.Frame(self)
        self.selectFrame.grid(row=1, column=0, sticky=tk.W+tk.E)

        # Tally selection
        labelTally = tk.Label(self.selectFrame, text='Tally:')
        labelTally.grid(row=0, column=0, sticky=tk.W)
        self.tallyBox = ttk.Combobox(self.selectFrame, state='readonly')
        self.tallyBox['values'] = [self.datafile.tallies[i].id
                                   for i in self.meshTallies]
        self.tallyBox.current(0)
        self.tallyBox.grid(row=0, column=1, sticky=tk.W+tk.E)
        self.tallyBox.bind('<<ComboboxSelected>>', self.update)

        # Planar basis selection
        labelBasis = tk.Label(self.selectFrame, text='Basis:')
        labelBasis.grid(row=1, column=0, sticky=tk.W)
        self.basisBox = ttk.Combobox(self.selectFrame, state='readonly')
        self.basisBox['values'] = ('xy', 'yz', 'xz')
        self.basisBox.current(0)
        self.basisBox.grid(row=1, column=1, sticky=tk.W+tk.E)
        self.basisBox.bind('<<ComboboxSelected>>', self.update)

        # Axial level
        labelAxial = tk.Label(self.selectFrame, text='Axial level:')
        labelAxial.grid(row=2, column=0, sticky=tk.W)
        self.axialBox = ttk.Combobox(self.selectFrame, state='readonly')
        self.axialBox.grid(row=2, column=1, sticky=tk.W+tk.E)
        self.axialBox.bind('<<ComboboxSelected>>', self.redraw)

        # Option for mean/uncertainty
        labelMean = tk.Label(self.selectFrame, text='Mean/Uncertainty:')
        labelMean.grid(row=3, column=0, sticky=tk.W)
        self.meanBox = ttk.Combobox(self.selectFrame, state='readonly')
        self.meanBox['values'] = ('Mean', 'Absolute uncertainty',
                                  'Relative uncertainty')
        self.meanBox.current(0)
        self.meanBox.grid(row=3, column=1, sticky=tk.W+tk.E)
        self.meanBox.bind('<<ComboboxSelected>>', self.update)

        # Scores
        labelScore = tk.Label(self.selectFrame, text='Score:')
        labelScore.grid(row=4, column=0, sticky=tk.W)
        self.scoreBox = ttk.Combobox(self.selectFrame, state='readonly')
        self.scoreBox.grid(row=4, column=1, sticky=tk.W+tk.E)
        self.scoreBox.bind('<<ComboboxSelected>>', self.redraw)

        # Filter label
        font = tkFont.Font(weight='bold')
        labelFilters = tk.Label(self.selectFrame, text='Filters:', font=font)
        labelFilters.grid(row=5, column=0, sticky=tk.W)

    def update(self, event=None):
        if not event:
            widget = None
        else:
            widget = event.widget

        tally_id = self.meshTallies[self.tallyBox.current()]
        selectedTally = self.datafile.tallies[tally_id]

        # Get mesh for selected tally
        self.mesh = self.datafile.meshes[
            selectedTally.filters['mesh'].bins[0] - 1]

        # Get mesh dimensions
        self.nx, self.ny, self.nz = self.mesh.dimension

        # Repopulate comboboxes baesd on current basis selection
        text = self.basisBox['values'][self.basisBox.current()]
        if text == 'xy':
            self.axialBox['values'] = [str(i+1) for i in range(self.nz)]
        elif text == 'yz':
            self.axialBox['values'] = [str(i+1) for i in range(self.nx)]
        else:
            self.axialBox['values'] = [str(i+1) for i in range(self.ny)]
        self.axialBox.current(0)

        # If update() was called by a change in the basis combobox, we don't
        # need to repopulate the filters
        if widget == self.basisBox:
            self.redraw()
            return

        # Update scores
        self.scoreBox['values'] = selectedTally.scores
        self.scoreBox.current(0)

        # Remove any filter labels/comboboxes that exist
        for row in range(6, self.selectFrame.grid_size()[1]):
            for w in self.selectFrame.grid_slaves(row=row):
                w.grid_forget()
                w.destroy()

        # create a label/combobox for each filter in selected tally
        count = 0
        for filterType in selectedTally.filters:
            if filterType == 'mesh':
                continue
            count += 1

            # Create label and combobox for this filter
            label = tk.Label(self.selectFrame, text=self.labels[filterType])
            label.grid(row=count+6, column=0, sticky=tk.W)
            combobox = ttk.Combobox(self.selectFrame, state='readonly')
            self.filterBoxes[filterType] = combobox

            # Set combobox items
            f = selectedTally.filters[filterType]
            if filterType in ['energyin', 'energyout']:
                combobox['values'] = ['{0} to {1}'.format(*f.bins[i:i+2])
                                      for i in range(f.length)]
            else:
                combobox['values'] = [str(i) for i in f.bins]

            combobox.current(0)
            combobox.grid(row=count+6, column=1, sticky=tk.W+tk.E)
            combobox.bind('<<ComboboxSelected>>', self.redraw)

        # If There are no filters, leave a 'None available' message
        if count == 0:
            count += 1
            label = tk.Label(self.selectFrame, text="None Available")
            label.grid(row=count+6, column=0, sticky=tk.W)

        self.redraw()

    def redraw(self, event=None):
        basis = self.basisBox.current() + 1
        axial_level = self.axialBox.current() + 1
        is_mean = self.meanBox.current()

        # Get selected tally
        tally_id = self.meshTallies[self.tallyBox.current()]
        selectedTally = self.datafile.tallies[tally_id]

        # Create spec_list
        spec_list = []
        for f in selectedTally.filters.values():
            if f.type == 'mesh':
                continue
            index = self.filterBoxes[f.type].current()
            spec_list.append((f.type, index))

        # Take is_mean and convert it to an index of the score
        score_loc = is_mean
        if score_loc > 1:
            score_loc = 1

        text = self.basisBox['values'][self.basisBox.current()]
        if text == 'xy':
            matrix = np.zeros((self.nx, self.ny))
            for i in range(self.nx):
                for j in range(self.ny):
                    matrix[i, j] = self.datafile.get_value(tally_id,
                        spec_list + [('mesh', (i + 1, j + 1, axial_level))],
                        self.scoreBox.current())[score_loc]
                    # Calculate relative uncertainty from absolute, if requested
                    if is_mean == 2:
                        # Take care to handle zero means when normalizing
                        mean_val = self.datafile.get_value(tally_id,
                            spec_list + [('mesh', (i + 1, j + 1, axial_level))],
                            self.scoreBox.current())[0]
                        if mean_val > 0.0:
                            matrix[i, j] = matrix[i, j] / mean_val
                        else:
                            matrix[i, j] = 0.0

        elif text == 'yz':
            matrix = np.zeros((self.ny, self.nz))
            for i in range(self.ny):
                for j in range(self.nz):
                    matrix[i, j] = self.datafile.get_value(tally_id,
                        spec_list + [('mesh', (axial_level, i + 1, j + 1))],
                        self.scoreBox.current())[score_loc]
                    # Calculate relative uncertainty from absolute, if requested
                    if is_mean == 2:
                        # Take care to handle zero means when normalizing
                        mean_val = self.datafile.get_value(tally_id,
                            spec_list + [('mesh', (axial_level, i + 1, j + 1))],
                            self.scoreBox.current())[0]
                        if mean_val > 0.0:
                            matrix[i, j] = matrix[i, j] / mean_val
                        else:
                            matrix[i, j] = 0.0

        else:
            matrix = np.zeros((self.nx, self.nz))
            for i in range(self.nx):
                for j in range(self.nz):
                    matrix[i, j] = self.datafile.get_value(tally_id,
                        spec_list + [('mesh', (i + 1, axial_level, j + 1))],
                        self.scoreBox.current())[score_loc]
                    # Calculate relative uncertainty from absolute, if requested
                    if is_mean == 2:
                        # Take care to handle zero means when normalizing
                        mean_val = self.datafile.get_value(tally_id,
                            spec_list + [('mesh', (i + 1, axial_level, j + 1))],
                            self.scoreBox.current())[0]
                        if mean_val > 0.0:
                            matrix[i, j] = matrix[i, j] / mean_val
                        else:
                            matrix[i, j] = 0.0

        # Clear the figure
        self.fig.clear()

        # Make figure, set up color bar
        self.axes = self.fig.add_subplot(111)
        cax = self.axes.imshow(matrix.transpose(), vmin=0.0, vmax=matrix.max(),
                               interpolation='none', origin='lower')
        self.fig.colorbar(cax)

        self.axes.set_xticks([])
        self.axes.set_yticks([])
        self.axes.set_aspect('equal')

        # Draw canvas
        self.canvas.draw()

    def get_file_data(self, filename):
        # Create StatePoint object and read in data
        self.datafile = StatePoint(filename)
        self.datafile.read_results()
        self.datafile.generate_stdev()

        # Find which tallies are mesh tallies
        self.meshTallies = []
        for itally, tally in enumerate(self.datafile.tallies):
            if 'mesh' in tally.filters:
                self.meshTallies.append(itally)

        if not self.meshTallies:
            tkMessageBox.showerror("Invalid StatePoint File",
                                   "File does not contain mesh tallies!")
            sys.exit(1)


if __name__ == '__main__':
    # Hide root window
    root = tk.Tk()
    root.withdraw()

    # If no filename given as command-line argument, open file dialog
    if len(sys.argv) < 2:
        filename = tkFileDialog.askopenfilename(title='Select statepoint file',
                                                initialdir='.')
    else:
        filename = sys.argv[1]

    if filename:
        # Check to make sure file exists
        if not os.path.isfile(filename):
            tkMessageBox.showerror("File not found",
                                   "Could not find regular file: " + filename)
            sys.exit(1)

        app = MeshPlotter(root, filename)
        root.deiconify()
        root.mainloop()
