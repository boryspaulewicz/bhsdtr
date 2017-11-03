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
bool fill_crit(int c, pen p = lightgray){
  fill(subpath(m, times(m, crit_x(c))[0], times(m, crit_x(c + 1))[0]) -- (crit_x(c + 1), 0) -- (crit_x(c), 0) -- cycle, p);
  return true;
}

fill_crit(1);
fill_crit(7);
draw(a, dashed+gray);
draw(b, dashed+gray);
int fs = 8;
for(int c = 1; c < nof_c; ++c){
  draw((crit_x(c),0) -- (crit_x(c), ypart(max(a)) + 0.5), dashed+gray);
  label(format("$c_%d$", c), (crit_x(c),0), S, fontsize(fs));
}
draw(m);
draw(evidence, Arrow(TeXHead), L = Label("$s$", EndPoint, fontsize(fs)));
label("mapping distribution $f$", position = (0, 0.5), fontsize(fs), filltype = Fill(white));
label("evidence distributions $p(s|stim)$", position = (0, ypart(max(b)) + 0.5), fontsize(fs), align = S, filltype = Fill(white));
// label("$p(s|stim=1)$", (-dprim/2, ypart(max(a)) * 0.6), fontsize(fs), filltype = Fill(white));
Label l = Label("$\gamma_2 =\ln\bigg(\frac{\int_{c_1}^{c_2} f(s) \,ds}{\int_{c_7}^{\infty} f(s) \, ds}\bigg)$",
      fontsize(fs));
// label(l, (crit_x(1.3), -1), align = E);
label(l, (crit_x(1.5), 2.1), filltype = Fill(white));
label("$\gamma_2$", (crit_x(1.5), -0.5), align = N, fontsize(fs), filltype = Fill(white));
label("$\gamma_8 = 0$", (crit_x(7.2), -0.5), align = NE, fontsize(fs), filltype = Fill(white));
