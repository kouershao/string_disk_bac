
void sort(double seq[3],int &i1,int &i2,int &i3) {
  if (seq[0] <= seq[1] && seq[1] <= seq[2]) {
    i1 = 0;
    i2 = 1;
    i3 = 2;
  }
  else if (seq[0] <= seq[2] && seq[2] <= seq[1]) {
    i1 = 0;
    i2 = 2;
    i3 = 1;
  }
  else if (seq[1] <= seq[0] && seq[0] <= seq[2]) {
    i1 = 1;
    i2 = 0;
    i3 = 2;
  }
  else if (seq[1] <= seq[2] && seq[2] <= seq[0]) {
    i1 = 1;
    i2 = 2;
    i3 = 0;
  }
  else if (seq[2] <= seq[0] && seq[0] <= seq[1]) {
    i1 = 2;
    i2 = 0;
    i3 = 1;
  }
  else if (seq[2] <= seq[1] && seq[1] <= seq[0]) {
    i1 = 2;
    i2 = 1;
    i3 = 0;
  }
}
