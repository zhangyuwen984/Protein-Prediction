/ATOMIC_POSITIONS / {
  do {
      getline
      if (NF>3){
      print $2, $3, $4
      }
  } while (NF > 3)
}
