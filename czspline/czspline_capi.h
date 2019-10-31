/*
 * C callable prototypes for czplsine subroutines
 */

#ifdef __cplusplus
extern "C"{
#endif

  void czspline_init1(int *n, int bcs[], int *ier);
  void czspline_set_axes1(int *n, double x[], int *ier);
  void czspline_set_ishermite1(int flag[2], int *ier);
  void czspline_set_bctypes1(int flag[2], int *ier);
  void czspline_set_bcvals1(double bcval1[], int *ier);
  void czspline_setup1(int *n, double f[], int *ier);
  void czspline_interp1(double *y, double *g, int *ier);
  void czspline_interp1_cloud(int *m, double *y, double *g, int *ier);
  void czspline_interp1_array(int *m, double *y, double *g, int *ier);
  void czspline_free1(int *ier);
  void czspline_save1(char const *fname, int *ier);
  void czspline_load1(char const *fname, int *ier);
  void czspline_isindomain1(double *y, int *ier);
  void czspline_isindomain1_cloud(int *m, double y[], int *ier);
  void czspline_isindomain1_array(int *m, double y[], int *ier);
  void czspline_gradient1(double *y, double *df, int *ier);
  void czspline_gradient1_cloud(int *m, double y[], double df[], int *ier);
  void czspline_gradient1_array(int *m, double y[], double df[], int *ier);
  void czspline_derivative1(int *i, double *y, double *f, int *ier);
  void czspline_derivative1_cloud(int *i, int *m, double y[], double f[], int *ier);
  void czspline_derivative1_array(int *i, int *m, double y[], double f[], int *ier);
  void czspline_isregular1(int *ier);

  void czspline_init2(int *n1, int *n2, int bcs1[], int bcs2[], int *ier);
  void czspline_set_axes2(int *n1, int *n2, double x1[], double x2[], int *ier);
  void czspline_set_ishermite2(int flag[], int *ier);
  void czspline_set_bctypes2(int flag[], int *ier);
  void czspline_set_bcvals2(double bcval1[], double bcval2[], int *ier);
  void czspline_setup2(int *n1, int *n2, double f[], int *ier);
  void czspline_interp2(double *y1, double *y2, double *g, int *ier);
  void czspline_interp2_cloud(int *m, double y1[], double y2[], double g[], int *ier);
  void czspline_interp2_array(int *m1, int *m2, double y1[], double y2[], double g[], int *ier);
  void czspline_free2(int *ier);
  void czspline_save2(char const *fname, int *ier);
  void czspline_load2(char const *fname, int *ier);
  void czspline_isindomain2(double *y1, double *y2, int *ier);
  void czspline_isindomain2_cloud(int *m, double y1[], double y2[], int *ier);
  void czspline_isindomain2_array(int *m1, int *m2, double y1[], double y2[], int *ier);
  void czspline_gradient2(double *y1, double *y2, double *df, int *ier);
  void czspline_gradient2_cloud(int *m, double y1[], double y2[], double df[], int *ier);
  void czspline_gradient2_array(int *m1, int *m2, double y1[], double y2[], double df[], int *ier);
  void czspline_derivative2(int *i1, int *i2, double *y1, double *y2, double *f, int *ier);
  void czspline_derivative2_cloud(int *i1, int *i2, int *m, double y1[], double y2[], double f[], int *ier);
  void czspline_derivative2_array(int *i1, int *i2, int *m1, int *m2, double y1[], double y2[], double f[], int *ier);
  void czspline_isregular2(int *ier);

  void czspline_init3(int *n1, int *n2, int *n3, int bcs1[], int bcs2[], int bc3[], int *ier);
  void czspline_set_axes3(int *n1, int *n2, int *n3, double x1[], double x2[], double x3[], int *ier);
  void czspline_set_ishermite3(int flag[], int *ier);
  void czspline_set_bctypes3(int flag[], int *ier);
  void czspline_set_bcvals3(double bcval1[], double bcval2[], double bcval3[], int *ier);
  void czspline_setup3(int *n1, int *n2, int *n3, double f[], int *ier);
  void czspline_interp3(double *y1, double *y2, double *y3, double *g, int *ier);
  void czspline_interp3_cloud(int *m, double y1[], double y2[], double y3[], double g[], int *ier);
  void czspline_interp3_array(int *m1, int *m2, int *m3, double y1[], double y2[], double y3[], double g[], int *ier);
  void czspline_free3(int *ier);
  void czspline_save3(char const *fname, int *ier);
  void czspline_load3(char const *fname, int *ier);
  void czspline_isindomain3(double *y1, double *y2, double *y3, int *ier);
  void czspline_isindomain3_cloud(int *m, double y1[], double y2[], double y3[], int *ier);
  void czspline_isindomain3_array(int *m1, int *m2, int *m3, double y1[], double y2[], double y3[], int *ier);
  void czspline_gradient3(double *y1, double *y2, double *df, int *ier);
  void czspline_gradient3_cloud(int *m, double y1[], double y2[], double y3[], double df[], int *ier);
  void czspline_gradient3_array(int *m1, int *m2, int *m3, double y1[], double y2[], double y3[], double df[], int *ier);
  void czspline_derivative3(int *i1, int *i2, int *i3, double *y1, double *y2, double y3[], double *f, int *ier);
  void czspline_derivative3_cloud(int *i1, int *i2, int *m, double y1[], double y2[], double y3[], double f[], int *ier);
  void czspline_derivative3_array(int *i1, int *i2, int *i3, int *m1, int *m2, int *m3, double y1[], double y2[], double y3[], double f[], int *ier);
  void czspline_isregular3(int *ier);

#ifdef __cplusplus
}
#endif
