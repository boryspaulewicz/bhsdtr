real normal(real x, real sigma)
{
  static real sqrt2pi=sqrt(2pi);
  return exp(-0.5*(x/sigma)^2)/(sigma*sqrt2pi);
}

settings.outformat = "pdf";

unitsize(1cm);
real dprim = 1.8;
real halfwidth = dprim + 2;
path evidence = (-(halfwidth),0)--(halfwidth + 0.1,0);
path a, b, m;
for(real x = -halfwidth; x <= halfwidth; x = x + 0.1){
  a = a -- (x, normal(x - dprim / 2, 1));
  b = b -- (x, normal(x + dprim / 2, 1));
  m = m -- (x, normal(x, 2));
}
real ys = 7;
a = yscale(ys) * a;
b = yscale(ys) * b;
m = yscale(ys) * m;

// Rysujemy kryteria
int nof_c = 8;
real crit_x(real c){
  if(c < nof_c){
    return (c - nof_c / 2) / nof_c * halfwidth * 2;
  }else{
    return xpart(max(m));
  }
}
draw(a, black);
draw(b, black);
int fs = 8;
for(int c = 1; c < nof_c; ++c){
  draw((crit_x(c),0) -- (crit_x(c), ypart(max(a)) + 0.5), dashed+gray);
  label(format("$c_%d$", c), (crit_x(c),0), S, fontsize(fs));
  label(format("$y=%d$", c), ((crit_x(c-1) + crit_x(c)) / 2, 2), S, fontsize(fs));
}
label(format("$y=%d$", 8), ((crit_x(8-1) + crit_x(8)) / 2, 2), S, fontsize(fs));
draw(evidence, Arrow(TeXHead), L = Label("$s$", EndPoint, fontsize(fs)));
label("evidence distributions $p(s|stim)$", position = (0, ypart(max(b)) + 0.5), fontsize(fs), align = S, filltype = Fill(white));