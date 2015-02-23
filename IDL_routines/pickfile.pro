; $Id: //depot/Release/ENVI52_IDL84/idl/idldir/lib/obsolete/pickfile.pro#1 $
;
; Copyright (c) 1991-2014, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.
;+
; NAME:
;       PICKFILE
;
; PURPOSE:
;       This function allows the user to interactively pick a file.  A file
;       selection tool with a graphical user interface is created.  Files
;       can be selected from the current directory or other directories.
;
; CATEGORY:
;       Widgets.
;
; CALLING SEQUENCE:
;       Result = PICKFILE()
;
; KEYWORD PARAMETERS:
;
;       FILE:   A string value for setting the initial value of the
;               selection. Useful if there is a default file
;
;       GET_PATH: Set to a named variable. Returns the path at the
;               time of selection.
;
;       GROUP:  The widget ID of the widget that calls PICKFILE.  When this
;               ID is specified, a death of the caller results in the death of
;               the PICKFILE widget application.
;
;       READ:   Set this keyword to make the title of the PICKFILE window
;               "Select File to Read".
;
;       WRITE:  Set this keyword to make the title of the PICKFILE window
;               "Select File to Write".
;
;       PATH:   The initial path to select files from.  If this keyword is
;               not set, the current directory is used.
;
;       FILTER: A string value for filtering the files in the file list.  This
;               keyword is used to reduce the number of files to choose from.
;               The user can modify the filter unless the FIX_FILTER keyword
;               is set.  Example filter values might be "*.pro" or "*.dat".
;
;       FIX_FILTER: When this keyword is set, only files that satisfy the
;               filter can be selected.  The user has no ability to modify
;               the filter and the filter is not shown.
;
;       TITLE:  A scalar string to be used for the window title.  If it is
;               not specified, the default title is "Select File"
;
;       NOCONFIRM: Return immediately upon selection of a file.  The default
;               behavior is to display the selection and then return the
;               file when the user uses the "ok" button.
;
;       MUST_EXIST: When set, only files that actually exist can be selected.
;
; OUTPUTS:
;       PICKFILE returns a string that contains the name of the file selected.
;       If no file is selected, PICKFILE returns a null string.
;
; COMMON BLOCKS:
;       NONE.
;
; SIDE EFFECTS:
;       NONE.
;
; RESTRICTIONS:
;       Only one instance of the PICKFILE widget can be running at one time.
;
; PROCEDURE:
;       Create and register the widget and then exit, returning the filename
;       that was picked.
;
; EXAMPLE:
;       Create a PICKFILE widget that lets users select only files with
;       the extensions 'pro' and 'dat'.  Use the 'Select File to Read' title
;       and store the name of the selected file in the variable F.  Enter:
;
;               F = PICKFILE(/READ, FILTER = '*.pro *.dat')
;
; MODIFICATION HISTORY:
;       Written by:     Steve Richards, April, 1991
;       July, 1991      Added a FILTER keyword to allow users
;                       to select files with a given extension or
;                       extensions.
;       August, 1991    Fixed bugs caused by differences between
;                       spawned ls commands on different machines.
;       September, 1991 Made Myfindfile so only one pass was
;                       necessary to find files and directories.
;       3/92 - ACY      Corrected initialization of dirsave, change spawn
;                       command to "ls -lL" and added case for links
;                       add NOCONFIRM keyword for auto exiting on selection
;       8/92 - SMR      Rewrote pickfile as a compound widget.
;       10/92 - SMR     Fixed a bug where extremely large file namess didn't
;                       show up properly in the file list or as return
;                       values.
;       12/92 - JWG     Add better machine dependency code
;       1/93 - JWG      Added FILE, GET_PATH keywords.
;       1/93 - TAC      Added Windows Common dialog pickfile code
;       2/93 - SMR      Fixed the documentation example for multiple extensions
;       1/94 - KDB      If directory had no execute permission on Unix
;                       platforms, CD fails and causes error. Added check
;                       for this. Increased spawn speed by using /sh for unix.
;                       Added -a switch to ls so that all files can be found
;                       on unix machines.
;       2/94 - KDB	Values passed to CD cannot end in a '\' on DOS
;			platforms. Program would crash if the PATH keyword
;			was supplied a value that ended with a "\". Added
;		        a check for this.
;	3/94 - BMH	Deleted the reference here to OS_PICKFILE for the
;			Unix platforms and created an IDL routine to
;			to call the Mac and Windows specific OS_PICKFILE
;			routines.  This solved the saving and restoring on
;	 		different platforms problem.
;	4/94 - KDB      The vms call to lib$findfile in valid_dir was
;		        commented out. This caused errors when path was
;			changed by user. Uncommented. In Valid_Dir, with
;			vms the type of directory specification was not
;		        checked (directory can be a path or a filename):
;			Fixed this. In dirlist section of event handler,
;		        a "[-]" would get trimmed to "" and cause error:
;			Fixed.
;	8/94 - ACY      Change the spawn command in getdirs to send error
;			output to /dev/null.
;	12/94 - DJE	Fix the FIX_FILTER option for the MacOS.
;	1/96 - RPM	Fixed reading of directories for when Unix long listing
;			(ls -l) does not align columns.
;	3/96 - LP	Implemented widget_pickfile in Motif
;                       (xmFileSelectionBox). Conformed to widget_...
;                       mechanism for all platforms.
;			Used Motif widget_pickfile, removed file, directory
;			internal IDL handling.
;	4/96 - LP	Renamed widget_pickfile to dialog_pickfile.
;
; Note. This routine is maintained for compatibility reason. New system
;       routine dialog_pickfile should be used instead.
;-
;

;------------------------------------------------------------------------------
;       procedure Pickfile
;------------------------------------------------------------------------------
FUNCTION Pickfile, Group=GROUP, Path=PATH, Read=READ, Write=WRITE, $
                   Filter=FILTER, Title=TITLE, Noconfirm=NOCONFIRM, $
                   Must_exist=MUST_EXIST, Fix_filter=FIX_FILTER, $
                   FILE=FILE, GET_PATH=GET_PATH
    
    here = ""
    thefile = ""
    filt = ""
    IF (KEYWORD_SET(NOCONFIRM))     THEN auto_exit = 1      ELSE auto_exit = 0
    IF (KEYWORD_SET(MUST_EXIST))    THEN existflag = 1      ELSE existflag = 0
    IF (KEYWORD_SET(FIX_FILTER))    THEN mapfilter = 1      ELSE mapfilter = 0
    IF (KEYWORD_SET(READ))          THEN readflag  = 1      ELSE readflag  = 0
    IF (KEYWORD_SET(WRITE))         THEN writeflag = 1      ELSE writeflag = 0
    IF (NOT (KEYWORD_SET(TITLE)))   THEN TITLE = ""
    IF (N_ELEMENTS(FILE) EQ 0)      THEN FILE = ""
    IF (N_ELEMENTS(GROUP) EQ 0)     THEN GROUP = 0
    
    CASE !Version.os OF
        'vms': begin
            separator       = ''
            IF (KEYWORD_SET(FILTER)) THEN filt = FILTER[0]
        end
; WINDOWS does NOT want a \ at the end of the directory
        'Win32': begin
            separator = ''
            IF (KEYWORD_SET(FILTER)) THEN filt = FILTER[0] ELSE filt = "*.*"
        end
        'MacOS': begin
            separator = ""
            IF (KEYWORD_SET(FILTER)) THEN filt = FILTER[0] ELSE filt = "*"
        end
        ELSE: begin
            separator       = '/'
            IF (KEYWORD_SET(FILTER)) THEN filt = FILTER[0]
            end
    ENDCASE
    
    CD, CURRENT = dirsave
    
    IF (N_ELEMENTS(PATH) EQ 0) THEN BEGIN
        
        PATH = dirsave + separator
        here = PATH
        
    ENDIF ELSE BEGIN
        
        ;; When on a Dos platform the argument to CD cannot end in a '\' unless
        ;; it is a root directory of a drive (ie C:\). Because of this, check
        ;; If we must remove the last character of PATH. -KDB 2/4/94
        
        IF((!Version.os EQ 'Win32')AND  $
           (Strpos(path, '\', Strlen(PATH)-1)NE -1))THEN  BEGIN
            IF(strlen(path) GT 3)THEN  $ ; Root dirs are 3 chars long.
              path = Strmid( path, 0, Strlen(path)-1)
        ENDIF
        
        IF(STRPOS(PATH, separator, STRLEN(PATH)- 1) EQ -1) AND $
          (PATH NE separator)THEN $
          PATH = PATH + separator
        
        CATCH, errorStatus
        if errorStatus NE 0 then begin
            here = ""
        endif else begin
            CD, PATH            ;if the user selected
            CATCH, /CANCEL
            here = PATH         ;a path then use it
        endelse
    ENDELSE
    
    thefile = DIALOG_PICKFILE(GROUP = GROUP, FILTER = filt, TITLE = TITLE, $
                              MUST_EXIST = existflag, FILE = FILE, $
                              NOCONFIRM = auto_exit, FIX_FILTER = mapfilter, $
                              GET_PATH = here, PATH = PATH, $
                              READ = readflag, WRITE = writeflag)
    
    CD, dirsave
    GET_PATH = here
    RETURN, thefile
END ;====================== end of Pickfile routine ===========================

