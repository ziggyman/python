procedure reduce_gtcmos(indirs)

string indirs     = "@dirs.list"    {prompt="list of directories to reduce"}
string *dirlist

begin

  file   infile

  # --- Erzeugen von temporaeren Filenamen
  infile = mktemp ("tmp")

# --- Umwandeln der Listen von Frames in temporaere Files
  if ((substr(indirs,1,1) == "@" && access(substr(indirs,2,strlen(indirs)))) || (substr(indirs,1,1) != "@" && access(indirs))){
    sections(indirs, option="root", > infile)
    dirlist = infile
  }
  else{
   if (substr(indirs,1,1) != "@"){
    print("reduce_gtcmos: ERROR: "//indirs//" not found!!!")
   }
   else{
    print("reduce_gtcmos: ERROR: "//substr(indirs,2,strlen(indirs))//" not found!!!")
   }
   return
  }

# --- build output filenames
  while (fscan (dirlist, in) != EOF){
    chdir(in)
    !ls > ls.list
  }

end
