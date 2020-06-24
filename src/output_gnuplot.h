




void ascii_output_3D(INT32 TL);
void write_ascii3D_field(CBlock* block, INT32 Field_Type, INT32 TL);

void write_field(CBlock* block, INT32 *lay, INT32 Field_Type, INT32 TL);

void gnuplot_write_inner(bool* Flag, INT32 xlay, INT32 ylay, INT32 zlay);
void gnuplot_write_grid(CBlock *temp_blk, INT32 xlay, INT32 ylay, INT32 zlay);

void gnuplot_xlayer(D_REAL* Field, char *field_name, INT32 c, INT32 xlay, INT32 step);
void gnuplot_ylayer(D_REAL* Field, char *field_name, INT32 c, INT32 ylay, INT32 step);
void gnuplot_zlayer(D_REAL* Field, char *field_name, INT32 c, INT32 zlay, INT32 step);
