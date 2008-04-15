/* Need a small wrapper for the fortran main program, since fortran
   code should not be used as the main program, see information
   elsewhere */
#ifdef __PGI
extern "C" int mcotac1d_();
int main()
{
    return mcotac1d_();
}
#else
extern "C" int mcotac1d();
int main()
{
    return mcotac1d();
}
#endif
