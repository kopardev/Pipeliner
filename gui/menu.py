
from tkinter import Menu
from tkinter.messagebox import showinfo

def build_in_progress( msg="" ) :
	if not msg :
		msg = 'Sorry! We are building this functionality'
	showinfo( 'Build In Progress:', msg )

def about():
    info="""
    CCBR Pipeliner
    Version 4.0.2
    April 29, 2020
    """
    showinfo("CCBR Pipeliner\nVersion 4.0.2",info)

def getHelp():
	info="""
	Tutorials and Documentation:
	https://github.com/CCBR/Pipeliner/wiki
	Email:
	CCBR_Pipeliner@mail.nih.gov
	"""
	showinfo("HELP",info)

def add_menubar( root ):
	menubar = Menu(root)

	#file menu
	filemenu = Menu(menubar, tearoff=0)
	# filemenu.add_command( label="Load project", command=build_in_progress )
	# filemenu.add_command( label="Save project", command=build_in_progress )
	filemenu.add_separator()
	filemenu.add_command( label="Exit", command=root.quit )

	menubar.add_cascade( label="File", menu=filemenu )

	#view menu
	viewmenu = Menu(menubar, tearoff=0)
	# viewmenu.add_command( label="Progress", command=root.progress )
	viewmenu.add_command( label="Workflow", command=root.workflow )

	menubar.add_cascade( label="View", menu=viewmenu )

	#help menu
	helpmenu = Menu(menubar, tearoff=0)
	helpmenu.add_command( label="GetHelp", command=getHelp )
	helpmenu.add_separator()
	helpmenu.add_command( label="About", command=about )

	menubar.add_cascade(label="Help", menu=helpmenu )
	
	root.config( menu=menubar )
