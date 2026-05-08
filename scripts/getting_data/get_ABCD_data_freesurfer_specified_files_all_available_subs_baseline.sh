# Script for downloading ABCD data freesurfer files (only desired files, for all available subjects) BASELINE - created by Bianca Serio on 21.10.24
# Modified on 10.06.25 to get the stats folders for Bin 
# Modified on 11.06.25: added a check to see if stats folder already was gotten (only proceed to download if empty) given that the script crashed in the middle of being run. Only made this modification for this script (not the equivalent for fu2y and fu4y)

# steps to run script: 
# (0. datalad clone ABCD data -> e.g., pt_02667/data/ABCD/ABCD_freesurfer)
# 1. check that directories defined in script match
# 2. ssh on MPI server -> getserver -sL
# 3. activate ssh agent to avoid needing to give pwd for datalad get each time ->  eval "$(ssh-agent -s)" then ssh-add ~/.ssh/id_ed25519_mpi 



#“#!” is an operator called shebang which directs the script to the interpreter location. So, if we use”#! /bin/sh” the script gets directed to the bourne-shell.
#!/bin/sh


# pptDIR is the participants directory
pptDIR=/data/pt_02667/data/ABCD/ABCD_freesurfer/

# defining the path to all ppt files 
pptFILES="/data/pt_02667/data/ABCD/ABCD_freesurfer/sub*"


# defining a counter for the downloaded subjects
count_subjects_downloaded=0



# LOOP OVER ALL PPT FILES
for pptfolder_path in $pptFILES
  do
    
    ### Given that there are some files with .html extensions, I want to filter those out (with isthereextension==0 -> no extension; 1 -> yes extension)
    ## obtain the number of fields in a path delimitated by .
    ## $filepath | tr -cd '/' | wc -c    -> tr -cd '.' removes all characters other than '.'; wc -c counts the remaining characters
    ## $(echo ) is used to provide an output to either expr command or to assign content to numfield variable
    isthereextension=$(echo $pptfolder_path | tr -cd '.' | wc -c)
    
    # if no extension (isthereextension==0)
    if [ "$isthereextension" -eq 0 ]; then
    
      ### To retrieve the subject name (instead of printing the whole path) -> truncate at every '/', count the number of fields in the whole path, +1 to take last one
      numfield=$(echo $(expr $(echo $pptfolder_path | tr -cd '/' | wc -c) + 1))
   
      ### Make a variable that contains just the file name - using the cut function 
      # -d (delimiter: using symbol "/" because path contains those)
      # -f (field: when you cut by fields (which are delimitated by a certain symbol, here .) you can chose the number of the field you want to select (defined by variable $numfield))
      ppt_name=$(echo $pptfolder_path | cut -d/ -f$numfield)
    
      
      
      echo "------------------Checking data for $ppt_name------------------"

      # Check if stats directory exists and is empty
      if [ -d "$pptfolder_path/stats" ] && [ -z "$(ls -A $pptfolder_path/stats)" ]; then

        echo "Stats folder exists and is empty. Proceeding with datalad get for $ppt_name..."

        cd $pptfolder_path

        # first do general datalad get --no-data (for all) which yields the broken symlinks in order for syntax with * to work to get specific files
        datalad get --no-data .
	
	# "datalad get X" for specific file structures
      	# # only use datalad get command once and concatenate the file name structures (will only log into JUDAC once and therefore be faster)
      	# # use -J NJOBS for faster execution (parallelizing via more nodes -> not sure if -J 10 is optimized)
      
     	#datalad get -J 10 surf/lh.area surf/rh.area  #previously
	
        datalad get -J 10 stats


      else
        echo "Skipping $ppt_name: either no stats folder or it is not empty."
      fi


      
      ### keeping track of number of subjects checked
      let "count_subjects_downloaded++"  # Increment the counter by 1
      echo "Number of subjects checked: $count_subjects_downloaded"
      
      
    fi

done

