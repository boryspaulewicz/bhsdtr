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

draw(a, black);
draw(b, black);
int fs = 8;
// Rysujemy pojedyncze kryterium
draw((0,0) -- (0, dist_top), dashed+gray);
label("$c$", (0, 0), S, fontsize(fs));
// Rysujemy etykiety rozk³adów
label("''noise''", position = (-dprim / 2, dist_top * .75), fontsize(fs * .5), align = S, filltype = Fill(white));
label("$(stim = 1)$", position = (-dprim / 2, dist_top * .65), fontsize(fs * .5), align = S, filltype = Fill(white));
label("''signal''", position = (dprim / 2, dist_top * .75), fontsize(fs * .5), align = S, filltype = Fill(white));
label("$(stim = 2)$", position = (dprim / 2, dist_top * .65), fontsize(fs * .5), align = S, filltype = Fill(white));
draw(evidence, Arrow(TeXHead), L = Label("$s$", EndPoint, fontsize(fs)));
// label("evidence distributions", position = (0, dist_top * .5), fontsize(fs), align = S, filltype = Fill(white));
draw((-dprim / 2, dist_top + .1) -- (dprim / 2, dist_top + .1), L = Label("$d'$", align = N, fontsize(fs)));
draw((-dprim / 2, dist_top + .1 + .05) -- (-dprim / 2, dist_top + .1 - .05));
draw((dprim / 2, dist_top + .1 + .05) -- (dprim / 2, dist_top + .1 - .05));
