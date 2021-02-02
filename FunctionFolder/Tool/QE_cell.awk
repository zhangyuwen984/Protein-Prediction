/CELL_PARAMETERS/ {
  for (i=1; i<=3; i++){
      getline
      print $1, $2, $3
  }
}
/crystal axes/ {
  for (i=1; i<=3; i++){
      getline
      print $4, $5, $6
  }
}
