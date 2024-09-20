library("bio3d")
generate_nma <- function(PDB, output){
  pdb <- read.pdb(PDB)
  modes <- nma(pdb)
  write(modes$fluctuations, file = output, sep="\n")
}
