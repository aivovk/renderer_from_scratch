template <class T, std::size_t n>
std::ostream& operator<<(std::ostream& out, const Vector<T,n>& r){
  out<<"("<<r.data[0];
  for(std::size_t i = 1 ; i < n ; i++)
    out<<","<<r.data[i];
  out<<")";
  return out;
}


template <class T, std::size_t nrows, std::size_t mcols>
std::ostream& operator<<(std::ostream& out, const Matrix<T, nrows, mcols>& m) {
  out<<"[";
  for(std::size_t row = 0 ; row < nrows ; row++) {
    if (row > 0)
      out << ' ';
    out << m.data[row][0];
    for(std::size_t col = 1 ; col < mcols ; col++)
      out<<" "<<m.data[row][col];
    if (row < nrows - 1)
      out<<'\n';
  }
  out<<"]";
  return out;

}
