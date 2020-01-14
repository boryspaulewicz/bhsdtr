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
real dist_top = ypart(max(a));

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
for(int c = 1; c < 3; ++c){
  draw((crit_x(c),0) -- (crit_x(c), dist_top), dashed+gray);
  label(format("$c_%d$", c), (crit_x(c),0), S, fontsize(fs));
  label(rotate(-90) * format("$y=%d$", c), ((crit_x(c-1) + crit_x(c)) / 2, 3), S, fontsize(fs));
}
draw((crit_x(nof_c - 1),0) -- (crit_x(nof_c - 1), dist_top), dashed+gray);
label("$c_{K-1}$", (crit_x(nof_c - 1),0), S, fontsize(fs));
label(rotate(-90) * "$y=K$", ((crit_x(nof_c - 1) + crit_x(nof_c)) / 2, 3), S, fontsize(fs));
draw((crit_x(nof_c - 2),0) -- (crit_x(nof_c - 2), dist_top), dashed+gray);
label("$c_{K-2}$", (crit_x(nof_c - 2),0), S, fontsize(fs));
label(rotate(-90) * "$y=K-1$", ((crit_x(nof_c - 2) + crit_x(nof_c - 1)) / 2, 3), S, fontsize(fs));
label("$\dots$", (crit_x(nof_c / 2),0), S, fontsize(fs));
// label("$\dots$", (crit_x(nof_c / 2), 1.9), S, fontsize(fs), filltype = Fill(white));
label("''noise''", position = (-dprim / 2, dist_top * .75), fontsize(fs * .5), align = S, filltype = Fill(white));
label("$(stim = 1)$", position = (-dprim / 2, dist_top * .65), fontsize(fs * .5), align = S, filltype = Fill(white));
label("''signal''", position = (dprim / 2, dist_top * .75), fontsize(fs * .5), align = S, filltype = Fill(white));
label("$(stim = 2)$", position = (dprim / 2, dist_top * .65), fontsize(fs * .5), align = S, filltype = Fill(white));
draw(evidence, Arrow(TeXHead), L = Label("$s$", EndPoint, fontsize(fs)));
// label("evidence distributions $p(s|stim)$", position = (0, ypart(max(b)) + 0.5), fontsize(fs), align = S, filltype = Fill(white));
