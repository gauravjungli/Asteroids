#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 16:25:58 2024

@author: g
"""

import tkinter as tk 
from tkinter import ttk
from PIL import Image, ImageTk  # For image loading (install Pillow if needed)
import sys
import os
from datetime import datetime
from Initialize import  Initialize_simulations
from IO import Output_File_old
import  subprocess
import threading
import queue
from plots import show_shape, show_omega, post_process, show_plots, show_3D_plots
from IO import Parameter
from tkinter import filedialog
from matplotlib.animation import  FFMpegWriter



class GUI:
    def __init__(self, screen_options):
        self.root = tk.Tk()
        self.root.protocol("WM_DELETE_WINDOW", self._quit)
        self.root.minsize(width=600, height=600)
        self.root.title("Asteroid simulation")
        self.root.withdraw()  # Hide the root window initially

        self.screen_options = screen_options
        self.parameters = {"run":0}
        self.Parameters()
        self.current_screen = 0
        self.screen_list =["welcome"] + list(self.screen_options.keys())+["preview","simulation","progress"]
        self.background_image = Image.open(os.path.join(os.path.dirname(os.path.dirname(os.getcwd())),"input","back_images.png"))  
        self.background_photo=None
        self.welcome_screen=None
        self.text_box = None
        self.widgets = None
        self.anim = None
        self.gui = None
        self.is_anim = False
        self.is_paused = False
        self.plot_index = -1
        self.next_frame = False
        self.show_next_screen()
        self.processes = {}  # Store subprocesses for termination
        #Uncomment this line only when you are trying to generate a parameters file for the solo run
        #Initialize_simulations(parameters=self.parameters)
    
    def Parameters(self):
        for name in self.screen_options:
            if name not in self.parameters:
                self.parameters[name]="Yes"
            inputs=self.screen_options[name]
            for Input in inputs:
                self.parameters[Input.Name]=Input.Value
        now=datetime.now()
        self.parameters['Output folder']= now.strftime("%Y-%m-%d_%H:%M")
    
    def _quit(self):
        self.root.quit()
        self.root.destroy() 
             
                           
    def show_doc(self, doc_text):
        if self.text_box:  
            self.text_box.delete("1.0", tk.END)
            self.text_box.insert(tk.END, doc_text)
    
    
    def create_doc_buttons(self, doc_list, frame):
    # Create the text box only once
        if self.text_box is None:
            self.text_box = tk.Text(self.root, wrap=tk.WORD, width=50, height=10)
            
            self.text_box.grid(row=0, column=4, rowspan=len(doc_list), sticky="nsew", padx=(5, 10), pady=5)
        for i, Input in enumerate(doc_list):
            button = tk.Button(frame, text="help", command=lambda text=Input.Help: self.show_doc(text))
            button.grid(row=i, column=3, sticky="w", padx=5, pady=5)
        
    def clear(self):
        for widget in self.root.winfo_children():
            widget.destroy()
        self.widgets={}

################################################################################

    def read_output(self,process, progress_queue):
        
        for line in iter(process.stdout.readline, ""):
    
            progress = float(line.strip())
            progress_queue.put(progress)
            print(f"Output from subprocess: {line}", end="")
                
    
    def run_script_instance(self,index,output_folder,run,progress_queue):
        python_path = sys.executable

        process = subprocess.Popen( [python_path, "main.py", 'Output folder', output_folder, 'run', str(run)], stdout=subprocess.PIPE,
        stderr  = subprocess.PIPE, text=True, bufsize=1 )
        
        #For debugging the code
        # Wait for the process to complete and get output/error
        # The timeout is optional but good practice
        #stdout_data, stderr_data = process.communicate(timeout=15)
        #return_code = process.returncode

        # print(f"Parent: Child STDOUT:\n{stdout_data}")
        #if stderr_data: # Only print if there's something in stderr
        #    print(f"Parent: Child STDERR:\n{stderr_data}")
        # print(f"Parent: Child Return Code: {return_code}")
        
        self.processes[index] = process
        self.read_output(process, progress_queue)


###############################################################################
        
    def update_progressbar(self, progress_bars, progress_queues, threads):
        for i, myqueue in enumerate(progress_queues):
            try:
                progress  = myqueue.get(block=False)
                progress_bars[i]["value"] = progress  # Update progress normally
            except:
                pass

                
    #schedule the next update after a short delay
        if all(not t.is_alive() for t in threads) or all(bar["value"] == 100 for bar in progress_bars):
            if  all(not t.is_alive() for t in threads):
                print("All threads ended")
            if all(bar["value"] == 100 for bar in progress_bars):
                print("All values in progress bar reached 100")
            text_frame = ttk.Frame(self.root)  # Create a custom frame style
            text_frame.pack(pady=(20, 10), padx=20)  # Adjust padding as needed
            
            # Software information
            ttk.Label(text_frame, text=r"Simulation Complete!", font=('Helvetica', 16, "bold")).pack(pady=10)
            ttk.Label(text_frame, text=r" Initializing post processing of the data", font=('Helvetica', 16, "bold")).pack(pady=10)
            if self.parameters['Dimension'] == '1D':
                post_process(self.parameters)
                ttk.Label(text_frame, text=" Post processing complete.", font=("Times", 16, "bold")).pack(pady=10)
                ttk.Label(text_frame, text=f"Outputs are written in  {self.parameters['Output folder']} folder ", font=("Times", 16, "bold")).pack(pady=10)
                if self.parameters["Single run"].lower() == 'no':
                    self.create_buttons(self.root,next_button_text="Show time evolution", back_button_text="",packing="pack")
            self.create_buttons(self.root,next_button_text="Start another simulation", back_button_text="Exit",packing="pack")

        else:      
            
            self.root.after(100, self.update_progressbar, progress_bars, progress_queues,threads)  
        
        
###################################################################
        
    def next_screen(self,text):
        
        
        if 0<self.current_screen<5:
            for key, widget in self.widgets.items():
                self.parameters[key] = widget.get()

        self.clear()
        
        if text == "Replay":
            self.current_screen -= 1
            
        if text == "Start another simulation":
            self.current_screen = 0
                    
        while True:

            self.current_screen+=1
            if  self.current_screen>4 or (self.parameters[self.screen_list[self.current_screen]]=="Yes" and self.parameters["Visualization"].lower()=="no"):
                break            
            
        self.show_next_screen()

#####################################################################

    def previous_screen(self,text):
        
        self.clear()
        
            
        if text == "Exit":
            self.current_screen=0
            

        while True:
            self.current_screen-=1
            if  self.current_screen<1 or self.parameters[self.screen_list[self.current_screen]]=="Yes":
                break            
        self.show_next_screen()
        
#######################################################################    
    
    def show_next_screen(self):

        
        if self.current_screen== -1:
            self.welcome_screen.destroy()
            self.root.destroy()
            print("Simulation aborted by the user")
            sys.exit()
            
        elif self.current_screen == 0:
            self.create_welcome_screen()
        
        elif self.current_screen < 5:
            if self.welcome_screen: 
                self.welcome_screen.destroy()
                self.root.deiconify()
            self.create_checklist_screen(self.screen_list[self.current_screen]) 
                
        elif self.current_screen == 5:
            
            if self.parameters["Visualization"].lower()=="yes":
                Parameter(self.parameters,"output")
                post_process(self.parameters) 
                self.create_post_processing_screen()
            else:
                self.create_preview_screen() 
        
        elif self.current_screen==6:
            self.create_simulation_screen()
            
        elif self.current_screen==7:
            self.create_progressbar_screen()

        elif self.current_screen==8:
            
            self.create_post_processing_screen()
            

        else:
            ttk.Label(self.root, text="No options selected.").pack(pady=20)

################################################################################
    
    def center_window(self,root):
        
        screen_width = root.winfo_screenwidth()
        screen_height = root.winfo_screenheight()
        width=root.winfo_reqwidth()
        height=root.winfo_reqheight()
        x = (screen_width ) // 4
        y = (screen_height ) // 3

        root.geometry(f"{width}x{height}+{x}+{y}")
        
        
################################################################################

    def disable_all_widgets(self,widget):
       
        if isinstance(widget, ttk.Entry) or isinstance(widget, ttk.Radiobutton):  # Check if it's a ttk widget
            widget.state(['disabled'])  
        for child in widget.winfo_children():  # Recursively disable children
            self.disable_all_widgets(child)
            
            
##################################################################################
    
    def create_buttons(self,root,packing="grid",next_button_text="Next",back_button_text="Back",row=1,column=0,
                       columnspan=1,myfont=('Helvetica', 16), *args):
        
        button_frame = ttk.Frame(root, height = 60)
        if self.current_screen==4:
            next_button_text="Preview"
        if packing=="grid":
            button_frame.grid(row=row, column=column,columnspan=columnspan, padx=10, pady=10, sticky="nsew")
        else:
            button_frame.pack( fill=tk.X)
        
        # Centering the buttons (same as before)
        button_frame.grid_columnconfigure(0, weight=1) 
        button_frame.grid_columnconfigure(1, weight=0) 
        button_frame.grid_columnconfigure(2, weight=0)
        button_frame.grid_columnconfigure(3, weight=1) 
        
        # Styling the buttons (corrected)
        style = ttk.Style()
        
        # Create new styles based on the default TButton style
        style.configure('Red.TButton', background='red', foreground='white',font=myfont)  # Use foreground for text color
        style.configure('Green.TButton', background='green', foreground='white', font=myfont)
            
        
        if back_button_text:
            back_button = ttk.Button(button_frame, text=back_button_text, command= lambda: self.previous_screen(back_button_text), style='Red.TButton')
            back_button.grid(row=0, column=1, padx=10, pady=10, sticky="ew")
        
        if next_button_text:
            next_button = ttk.Button(button_frame, text=next_button_text, command=lambda: self.next_screen(next_button_text), style='Green.TButton')
            next_button.grid(row=0, column=2, padx=10, pady=10, sticky="ew")
        
        
#######################################################################################
        

    def display_buttons(self,myfont=('Helvetica', 16)):
        
        button_frame = ttk.Frame(self.root)
        button_frame.pack(side=tk.BOTTOM)
        
        style = ttk.Style()
        
        # Create new styles based on the default TButton style
        style.configure('Red.TButton', background='red', foreground='white',font=myfont)  # Use foreground for text color
        style.configure('Green.TButton', background='green', foreground='white', font=myfont)
        
        self.play_pause_button = ttk.Button(button_frame, text="Pause", command=self.toggle_play_pause,style='Red.TButton')
        self.play_pause_button.pack(side=tk.LEFT,padx=5)
        
        self.restart_button = ttk.Button(button_frame, text="Restart", command=self.restart,style='Green.TButton')
        self.restart_button.pack(side=tk.LEFT,padx=5)
        
        self.previous_button = ttk.Button(button_frame, text="Previous", command=self.backward,style='Green.TButton')
        self.previous_button.pack(side=tk.LEFT)
        
        self.next_button = ttk.Button(button_frame, text="Next", command=self.forward,style='Green.TButton')
        self.next_button.pack(side=tk.LEFT)
        
        self.save_button = ttk.Button(button_frame, text="Save", command=self.save_animation,style='Green.TButton')
        self.save_button.pack(side=tk.LEFT,padx=5)
        
        
        self.restart_button = ttk.Button(button_frame, text="Start new simulation", command=lambda: self.next_screen("Start another simulation"),style='Green.TButton')
        self.restart_button.pack(side=tk.LEFT,padx=5)
        
        self.exit_button = ttk.Button(button_frame, text="Exit", command=self.root.destroy,style='Red.TButton')
        self.exit_button.pack(side=tk.LEFT,padx=10)
        
        
        
    def toggle_play_pause(self):
        if self.is_paused:
            self.is_paused = False
            self.play_pause_button.config(text="Pause",style='Red.TButton')
        else:
            self.is_paused = True
            self.play_pause_button.config(text="Play",style='Green.TButton')
            
        
    def restart(self):
        
        self.is_paused = False
        self.plot_index = -1
        self.play_pause_button.config(text="Pause",style='Red.TButton')

    def forward(self):
        
        self.next_frame = True
        self.plot_index += 1 
        
            
    def backward(self):
        
        self.next_frame = True
        self.plot_index -= 1
        
    def save_animation(self):
        
        self.plot_index = -1
        self.is_paused = True
        self.play_pause_button.config(text="Play",style='Green.TButton')
        self.is_anim = True
        self.progress_window = tk.Toplevel(self.root)
        self.progress_window.title("Saving Animation")
        
        self.progress_label = tk.Label(self.progress_window, text="Saving animation, please wait...")
        self.progress_label.pack(pady=10)
        
        self.progress = ttk.Progressbar(self.progress_window, orient='horizontal', length=200, mode='determinate')
        self.progress.pack(pady=10)
        self.progress_window.geometry("300x100+1000+500")
        filename=Output_File_old(self.parameters,"output")
        file_path = filedialog.asksaveasfilename ( initialdir =filename,
           defaultextension = ".mp4", filetypes = [("MP4 files", "*.mp4")],
           initialfile = "animation.mp4" )
        
        if file_path:
            writer = FFMpegWriter(fps=30, metadata=dict(artist='Me'), bitrate=15000)
            with writer.saving(self.fig, file_path, dpi=300 ):
                for i in range(self.progress['maximum'] ):
                    self.anim._draw_frame(i)
                    writer.grab_frame()
                    self.progress['value'] = i + 1
                    self.progress_window.update_idletasks()
            
            self.progress_window.destroy()
            print("Animation saved as animation.mp4!")


            self.anim.event_source.stop()
            self.anim = None

            

#########################################################################################

    def load_frame(self,frame,inputs):
        
        i=0
        for Input in inputs:
            
            label = ttk.Label(frame, text=Input.Name)
            label.grid(row=i, column=0, sticky="w")
            
            var=None
            if Input.Type=="bool":
    
             # ttk Styling (modified)
                style = ttk.Style()
                
                # Increase the indicator size slightly
                style.configure("TRadiobutton", indicatorsize=20) 
                
                # Change the indicator color when selected to create a circle effect
                style.map("TRadiobutton")   
                var = tk.StringVar(value=self.parameters[Input.Name])
                yes_radio = ttk.Radiobutton(frame, text="Yes", variable=var, value="Yes")
                no_radio = ttk.Radiobutton(frame, text="No", variable=var, value="No")
    
                # Place the radio buttons in the same grid location as the checkbox
                yes_radio.grid(row=i, column=1,  sticky="w")  # Adjust row/column as needed
                no_radio.grid(row=i,  column=2, sticky="w")

            elif Input.Type=="str" and Input.Options:
                var = tk.StringVar(value=self.parameters[Input.Name])
                combobox = ttk.Combobox(frame, values=Input.Options, textvariable=var)
                combobox.grid(row=i, column=1, columnspan=2, rowspan=1, sticky="nsew", padx=5, pady=5)
                #combobox.set(Input.Value)
            else:
    
                var = ttk.Entry(frame)
                var.insert(0,self.parameters[Input.Name])
                var.grid(row=i, column=1, columnspan=2, padx=5, pady=5)
            self.widgets[Input.Name]=var
            i+=1
            
##############################################################################################
    
    def create_checklist_screen(self,name):  
        inputs=self.screen_options[name]
    
        frame = ttk.LabelFrame(self.root, text=name)
        frame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
        # Text box for documentation
        
        self.text_box = tk.Text(self.root, wrap=tk.WORD, width=50, height=10)
        self.text_box.grid(row=0, column=4, rowspan=len(inputs), sticky="nsew", padx=(5, 10), pady=5)
        
        # Create the buttons
        self.create_doc_buttons(inputs, frame)
    
        # Configure grid layout
        self.root.columnconfigure(4, weight=1)
        self.root.rowconfigure(list(range(len(inputs))), weight=1)  # Expand rows
        
        self.load_frame(frame, inputs)    
        self.create_buttons(self.root,next_button_text="Next",back_button_text="Back")
        #self.center_window(self.root)

##############################################################################################    

    def create_simulation_screen(self):
        ttk.Label(self.root, text="Simulation running...").pack(pady=20)
        #self.center_window(self.root)
        self.root.after(1000, self.next_screen,"next")  # Simulate 2 seconds of running
        
#############################################################################################       
    
    def create_preview_screen(self):

              # --- 1. Parent Frame & Master Grid Weights ---
        main_frame = ttk.Frame(self.root)
        
        # Use pack to force the main_frame to fill the entire application window
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Keep the grid weights for the inside of main_frame so the Canvas expands
        main_frame.grid_rowconfigure(0, weight=1)    
        main_frame.grid_columnconfigure(0, weight=1)
        
        # --- 2. Scrollbars and Canvas ---
        x_scrollbar = ttk.Scrollbar(main_frame, orient=tk.HORIZONTAL)
        y_scrollbar = ttk.Scrollbar(main_frame, orient=tk.VERTICAL)
        
        my_canvas = tk.Canvas(main_frame, xscrollcommand=x_scrollbar.set, yscrollcommand=y_scrollbar.set, highlightthickness=0)
        
        my_canvas.grid(row=0, column=0, sticky="nsew")
        y_scrollbar.grid(row=0, column=1, sticky="ns")
        x_scrollbar.grid(row=1, column=0, sticky="ew")
        
        y_scrollbar.config(command=my_canvas.yview)
        x_scrollbar.config(command=my_canvas.xview)
        
        # --- 3. The Inner Scrollable Frame ---
        scrollable_frame = ttk.Frame(my_canvas)
        
        # CRITICAL FIX: Tell the inner frame's columns to expand to fill space!
        scrollable_frame.grid_columnconfigure(0, weight=1)
        scrollable_frame.grid_columnconfigure(1, weight=1)
        
        # Create window inside canvas
        canvas_window = my_canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        
        # Update scrollregion when inner frame changes size
        scrollable_frame.bind(
            "<Configure>",
            lambda e: my_canvas.configure(scrollregion=my_canvas.bbox("all"))
        )
        
        # Force the inner frame to match the Canvas dimensions
        def _on_canvas_configure(event):
            # Match Canvas width
            my_canvas.itemconfig(canvas_window, width=event.width)
            
            # Match Canvas height if content is small, otherwise allow scrolling (height=0)
            if scrollable_frame.winfo_reqheight() < event.height:
                my_canvas.itemconfig(canvas_window, height=event.height)
            else:
                my_canvas.itemconfig(canvas_window, height=0)
        
        my_canvas.bind("<Configure>", _on_canvas_configure)
        
        # --- 4. Your Dynamic Grid Loop ---
        j = 0
        column = 0
        
        for name in self.screen_options:
            if self.parameters[name] == "Yes":
                inputs = self.screen_options[name]
                
                # Parent is the inner scrollable_frame
                frame = ttk.LabelFrame(scrollable_frame, text=name)
                frame.grid(row=j, column=column, padx=10, pady=10, sticky="nsew")
                
                # Make the inside of the LabelFrame expand too (if needed for your widgets)
                frame.grid_columnconfigure(0, weight=1)
                frame.grid_columnconfigure(1, weight=1)
                
                # Alternating logic
                j += column
                column = int(not bool(column))
        
                # Load content into the frame
                self.load_frame(frame, inputs)
                self.disable_all_widgets(frame)
        
        # --- 5. Mouse Wheel Scrolling (Linux/Cross-platform) ---
        def _on_arrow_keys(event):
            if event.keysym == 'Up':
                my_canvas.yview_scroll(-1, "units")
            elif event.keysym == 'Down':
                my_canvas.yview_scroll(1, "units")
        
        my_canvas.bind_all("<Button-4>", lambda e: _on_arrow_keys(e) or my_canvas.yview_scroll(-1, "units")) # Linux Up
        my_canvas.bind_all("<Button-5>", lambda e: _on_arrow_keys(e) or my_canvas.yview_scroll(1, "units")) # Linux Down


        self.create_buttons(self.root,next_button_text="Start simulation", back_button_text="Back",packing="pack",row=1,column=0,columnspan=2)
#self.center_window(self.root)
 
###############################################################################################        

            
    def create_progressbar_screen(self):

        # --- Parent Frame for Everything ---
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=1)
        
        # --- NEW: A frame to hold the canvas and the VERTICAL scrollbar together ---
        top_container = ttk.Frame(main_frame)
        # This container will expand to fill the space *above* the horizontal scrollbar
        top_container.pack(fill=tk.BOTH, expand=1)
        
        # --- The Horizontal scrollbar goes at the very bottom of the main_frame ---
        x_scrollbar = ttk.Scrollbar(main_frame, orient=tk.HORIZONTAL)
        x_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)
        
        # --- The Canvas and Vertical scrollbar go INSIDE the top_container ---
        my_canvas = tk.Canvas(top_container, xscrollcommand=x_scrollbar.set)
        my_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        
        y_scrollbar = ttk.Scrollbar(top_container, orient=tk.VERTICAL, command=my_canvas.yview)
        y_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # --- Final Configuration ---
        my_canvas.configure(yscrollcommand=y_scrollbar.set)
        x_scrollbar.config(command=my_canvas.xview) # Configure x-scrollbar command here
        
        # Bind and add the inner frame
        my_canvas.bind('<Configure>', lambda e: my_canvas.configure(scrollregion = my_canvas.bbox("all")))
        scrollable_frame = ttk.Frame(my_canvas)
        my_canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        
        def _on_arrow_keys(event):
            # For Linux arrow key scrolling
            if event.keysym == 'Up':
                my_canvas.yview_scroll(-1, "units")
            elif event.keysym == 'Down':
                my_canvas.yview_scroll(1, "units")
        
        my_canvas.bind_all("<Button-4>", lambda e: _on_arrow_keys(e) or my_canvas.yview_scroll(-1, "units")) # Linux
        my_canvas.bind_all("<Button-5>", lambda e: _on_arrow_keys(e) or my_canvas.yview_scroll(1, "units")) # Linux
        
        parameters_list=[]
        Initialize_simulations(parameters=self.parameters,parameters_list=parameters_list)
        
        threads = []
        progress_bars = []
        progress_queues = [queue.Queue() for _ in parameters_list]
        for i, par in enumerate(parameters_list):
            frame = ttk.Frame(scrollable_frame)
            frame.pack(fill="x")  # Expand horizontally
            style = ttk.Style()
            style.configure("Bold.TLabel", font=("Helvetica", 12, "bold")) 
            label = ttk.Label(frame, text=f"run {par['run']}",  style="Bold.TLabel")
            label.pack(side="left",padx=10)
            progress_bar =  ttk.Progressbar(frame, orient="horizontal", length=300, mode="determinate")
            progress_bar.pack(side="left",padx=5,pady=10)
            # Individual abort button
            abort_button = ttk.Button(frame, text="Abort", command=lambda i=i: self.abort_process(i))
            abort_button.pack(side="right", padx=5)
            show_plot_button = ttk.Button(frame, text="Show shape", command=lambda par=par: show_shape(par))
            show_plot_button.pack(side="right", padx=5) 
            show_plot_button = ttk.Button(frame, text="Show spin", command=lambda par=par: show_omega(par))
            show_plot_button.pack(side="right", padx=5) 
            
            progress_bars.append(progress_bar)
            # Create an abort flag for the process

            thread = threading.Thread(target=self.run_script_instance, args=(i,par['Output folder'],par['run'],progress_queues[i]))
            thread.start()
            threads.append(thread)
        
        # Global abort button
        abort_all_button = ttk.Button(self.root, text="Abort Simulation", command=self.abort_all_processes)
        abort_all_button.pack(pady=10)
        self.update_progressbar(progress_bars, progress_queues, threads)
        
    def abort_process(self, index):
        """Abort a specific process"""
  
        if  self.processes[index].poll() is None:
            self.processes[index].terminate()
            print(f"Process {index + 1} terminated.")

    def abort_all_processes(self):
        """Abort all processes"""

        for i in self.processes:
            if self.processes[i].poll() is None:
                self.processes[i].terminate()
                print(f"Process {i + 1} terminated.")
        
######################################################################################################

    def create_post_processing_screen(self):
        
        self.display_buttons() 
        show_3D_plots(self)
        
    
##############################################################################################################
    """ Creates welcome screen"""
    
    def create_welcome_screen(self):
        
        self.welcome_screen = tk.Toplevel(self.root)
        self.welcome_screen.protocol("WM_DELETE_WINDOW", self._quit)
        self.welcome_screen.title("Welcome")
        self.welcome_screen.minsize(width=600, height=600)  # Adjust size as needed

        #Add welcome screen content (labels, images, etc.)
        welcome_label = ttk.Label(self.welcome_screen, text="Welcome to the Asteroid Simulation!")
        welcome_label.pack(pady=10)
        
        # Load and resize background image
        
        
        
        self.background_photo = ImageTk.PhotoImage(self.background_image)
        background_label = ttk.Label(self.welcome_screen, image=self.background_photo)
        background_label.place(x=0, y=0, relwidth=1, relheight=1)  # Cover the entire window
        
        
        #Load and display image (replace 'your_image.png' with your actual image path)
        image = Image.open(os.path.join(os.path.dirname(os.path.dirname(os.getcwd())),"input","images.png")) 
        image.thumbnail((900, 450)) 
        photo = ImageTk.PhotoImage(image)
        
        self.welcome_screen.image = photo
        image_label = ttk.Label(self.welcome_screen, image=photo)
        image_label.pack(pady=10)
        
        
        text_frame = ttk.Frame(self.welcome_screen)  # Create a custom frame style
        text_frame.pack(pady=(20, 10), padx=20)  # Adjust padding as needed
        
        # Software information
        ttk.Label(text_frame, text=r"SPASE", font=("Times", 46, "bold")).pack()
        ttk.Label(text_frame, text=r" Stochastic Program for Surface Evolution of Asteroids", font=("Times", 28, "bold")).pack()
        # ... (creator information remains the same, but place everything in text_frame instead of welcome_screen)
        
    
        # #Next button to proceed to checklists
        self.create_buttons(self.welcome_screen,next_button_text="Start", back_button_text="Abort",packing="pack", myfont=('Helvetica', 24))
        
        # Creators frame
        ttk.Label(self.welcome_screen, text="Created by:", font=("Helvetica", 12)).pack(pady=5)
        creators_frame = ttk.Frame(self.welcome_screen, style="Transparent.TLabel")
        creators_frame.pack()
        
        creators = [
            ("Gaurav", os.path.join(os.path.dirname(os.path.dirname(os.getcwd())),"input", "gaurav.jpeg")),
            ("Deepayan", os.path.join(os.path.dirname(os.path.dirname(os.getcwd())),"input", "Deepayan.jpeg")),
            ("Ishan", os.path.join(os.path.dirname(os.path.dirname(os.getcwd())),"input", "Ishan.jpeg")),
        ]
        
        max_width = 100
        max_height = 100
        
        for name, photo_path in creators:
            creator_subframe = ttk.Frame(creators_frame)
            creator_subframe.pack(side=tk.LEFT, padx=10)
        
            # Load image and maintain aspect ratio
            img = Image.open(photo_path)
            img.thumbnail((max_width, max_height))  # Resize while keeping aspect ratio
        
            creator_photo = ImageTk.PhotoImage(img)  # Create the PhotoImage from the resized PIL image
        
            photo_label = ttk.Label(creator_subframe, image=creator_photo)
            photo_label.image = creator_photo  # Keep a reference
            photo_label.pack()
        
            ttk.Label(creator_subframe, text=name).pack()
        
        #self.center_window(self.welcome_screen)

