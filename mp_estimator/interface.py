# -*- coding: utf-8 -*-
"""
    interface.py contains all functions to create the user interface to load
    the experimental data (that is: show each QR, SR and ISOM experiment and
    show the data that is selected)

    author: Edwin D.H.M. Reuvers
    v1.0.0, september 2025
"""

#%% Import packages
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

#%% Function to let user select data
def OpenFileDialog():
    # Create a Tkinter root window (it will not be displayed)
    root = tk.Tk()
    root.withdraw()  # Hide the root window

    # Open the file dialog to select multiple files
    file_paths = filedialog.askopenfilenames(
        title="Select Files",          # Title of the dialog
        filetypes=[("All Files", "*.*")]  # Filter for file types (optional)
    )
    
    return file_paths

#%% Function to select data column
def SelDataCol(nCol,nColMax):
    # Create a Tkinter root window
    root = tk.Tk()
    root.title("Select Numbers")

    # Set the size of the window
    root.geometry("500x400")  # Adjusted Width x Height

    # Variables to store the answers
    cols = {}
    is_submitted = [False]  # List to use as a mutable flag

    # Define a function to handle the submission
    def on_submit():
        try:
            selected_values = [int(combo.get()) for combo in comboboxes]
        except ValueError:
            messagebox.showerror("Input Error", "Please select valid numbers.")
            return

        # Check if all answers are within the valid range
        if any(val not in range(nColMax) for val in selected_values):
            messagebox.showerror("Input Error", "Please select numbers within the valid range.")
            return

        # Check for uniqueness of the answers
        if len(selected_values) != len(set(selected_values)):
            messagebox.showerror("Input Error", "Each number must be unique.")
            return

        # Update the answers dictionary
        for i, val in enumerate(selected_values, 1):
            cols[f'col{i}'] = val

        # Set the flag to indicate submission
        is_submitted[0] = True

        # Close the Tkinter window
        root.quit()
        root.destroy()

    # Define a function to handle the window close event
    def on_close():
        if not is_submitted[0]:  # Check if not submitted
            messagebox.showerror("Input Error", "Please submit the form before closing.")
        else:
            root.quit()
            root.destroy()

    # Bind the close button (X) of the window to the on_close function
    root.protocol("WM_DELETE_WINDOW", on_close)

    # Configure grid column and row weights
    for i in range(nCol + 1):
        root.grid_rowconfigure(i, weight=1)
    root.grid_columnconfigure(0, weight=1)
    root.grid_columnconfigure(1, weight=1)

    # Create and place labels and comboboxes dynamically based on nCol
    comboboxes = []
    labels = ['time', 'lmtc', 'fsee', 'stim']  # Extend as needed for more variables
    for i in range(nCol):
        tk.Label(root, text=f"Select the column number of variable: '{labels[i]}':", wraplength=250).grid(row=i, column=0, padx=20, pady=20, sticky="ew")
        combo = ttk.Combobox(root, values=list(range(nColMax)), state="readonly")
        combo.grid(row=i, column=1, padx=10, pady=10, sticky="ew")
        combo.set("Column number")
        comboboxes.append(combo)

    # Create a submit button
    submit_button = tk.Button(root, text="Submit", command=on_submit)
    submit_button.grid(row=nCol, column=0, columnspan=2, pady=10, sticky="ew")

    # Start the Tkinter event loop
    root.mainloop()

    # Return the answers after the Tkinter window is closed
    if not is_submitted[0]:
        raise RuntimeError("The form was not submitted before closing.")
    return tuple(cols[f'col{i+1}'] for i in range(nCol))

#%% Function to plot selected data
def DataPlot(time,lmtc,fsee,initial_ranges,filename=""):
    # Create a root Tkinter window (hidden)
    root = tk.Tk()
    root.withdraw()  # Hide the root window

    # Get screen dimensions
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()

    # Calculate sizes and positions
    plot_window_width = int(screen_width * 0.8)
    plot_window_height = int(screen_height * 0.8)
    feedback_window_width = screen_width - plot_window_width
    feedback_window_height = plot_window_height

    # Create a new Tkinter popup window for the plot
    popup_window = tk.Toplevel()
    popup_window.title("Plot Popup")
    popup_window.geometry(f"{plot_window_width}x{plot_window_height}+0+0")
    popup_window.attributes('-topmost', True)  # Keep the window on top

    # Create a new Tkinter popup window for feedback
    feedback_window = tk.Toplevel()
    feedback_window.title("Provide Feedback")
    feedback_window.geometry(f"{feedback_window_width}x{feedback_window_height}+{plot_window_width}+0")
    feedback_window.attributes('-topmost', True)  # Keep the window on top

    # Create a matplotlib figure and axis
    fig = Figure()
    fig.suptitle('File = '+ filename)
    ax1 = fig.add_subplot(211) 
    ax2 = fig.add_subplot(212, sharex=ax1)
    
    # Create a canvas to embed the plot in Tkinter
    canvas = FigureCanvasTkAgg(fig, master=popup_window)
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # Add the Matplotlib toolbar for zoom and pan
    toolbar = NavigationToolbar2Tk(canvas, popup_window)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    # Define current ranges and function to update plot
    current_ranges = initial_ranges.copy()

    def update_plot():
        axs = [ax1, ax2]
        for iAx, y in enumerate([lmtc,fsee]):
            ax = axs[iAx]
            # Clear the figure on top and plot data
            ax.clear()
            ax.plot(time, y)
    
            # Plot data based on updated ranges
            for param, (start, end) in current_ranges.items():
                start_idx, end_idx = int(start) - 1, int(end)
                if 0 <= start_idx < len(time) and 0 <= end_idx <= len(time):
                    ax.plot(time[start_idx:end_idx], y[start_idx:end_idx], '.', label=f'{param}')
    
            ax.set_xlabel('Time [s]')
            if iAx == 0:
                ax.set_ylabel('$L_{MTC}$ [m]')
            elif iAx == 1:
                ax.set_ylabel('$F_{SEE}$ [N]')
            ax.legend()
            canvas.draw()
               
    def on_close():
        # Output the current ranges before closing the window
        print("Figure closed by user, so current range is used. Which is:")
        for param, (start, end) in current_ranges.items():
            print(f"{param}: start={start}, end={end}")

        # Close the feedback window if it's open
        if feedback_window:
            feedback_window.destroy()

        popup_window.destroy()  # Close the popup window
        root.quit()  # End the Tkinter main loop

    # Bind the close event to our handler
    popup_window.protocol("WM_DELETE_WINDOW", on_close)
    feedback_window.protocol("WM_DELETE_WINDOW", on_close)

    # Define the ranges_entries at the outer scope
    ranges_entries = {}

    def preview_changes():
        # Update the ranges from the feedback window
        for param, (start_entry, end_entry) in ranges_entries.items():
            try:
                start = int(start_entry.get())
                end = int(end_entry.get())
                if start < 1 or end > len(time) or start > end:
                    raise ValueError
                current_ranges[param] = (start, end)
            except ValueError:
                messagebox.showerror("Error", f"Invalid range for {param}. Please enter valid start and end indices.")
                return

        update_plot()  # Update plot after getting new ranges

    def confirm_changes():
        preview_changes()  # Ensure the latest ranges are saved

        # Close both windows
        feedback_window.destroy()
        popup_window.destroy()
        root.quit()  # End the Tkinter main loop

    def provide_feedback():
        tk.Label(feedback_window, text="Set ranges for all parameters:").pack(pady=10)

        parameters = current_ranges.keys()
        for param in parameters:
            frame = tk.Frame(feedback_window)
            frame.pack(pady=5, fill=tk.X)
            
            tk.Label(frame, text=param, width=20, anchor='w').pack(side=tk.LEFT, padx=10)
            
            tk.Label(frame, text="Start:").pack(side=tk.LEFT)
            start_entry = tk.Entry(frame, width=5)
            start_entry.pack(side=tk.LEFT, padx=5)
            
            tk.Label(frame, text="End:").pack(side=tk.LEFT)
            end_entry = tk.Entry(frame, width=5)
            end_entry.pack(side=tk.LEFT, padx=5)
            
            # Set initial values
            start_value, end_value = current_ranges.get(param, (1, len(time)))
            start_entry.insert(0, start_value)
            end_entry.insert(0, end_value)
            
            ranges_entries[param] = (start_entry, end_entry)
        
        # Add buttons in a separate frame
        button_frame = tk.Frame(feedback_window)
        button_frame.pack(fill=tk.X, padx=10, pady=10)
        
        preview_button = tk.Button(button_frame, text="Preview changes", command=preview_changes)
        preview_button.pack(side=tk.LEFT, padx=5)
    
        confirm_button = tk.Button(button_frame, text="Confirm changes", command=confirm_changes)
        confirm_button.pack(side=tk.LEFT, padx=5)

    # Show initial plot and feedback
    update_plot()
    provide_feedback()

    # Start the Tkinter event loop
    root.mainloop()

    # Return both initial and final ranges
    return current_ranges

