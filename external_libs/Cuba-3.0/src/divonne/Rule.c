/*
	Rule.c
		integration with cubature rules
		code lifted with minor modifications from DCUHRE
		by J. Berntsen, T. Espelid, and A. Genz
		this file is part of Divonne
		last modified 17 Dec 11 th
*/


enum { nrules = 5 };

#define TYPEDEFSET \
  typedef struct { \
    count n; \
    real weight[nrules], scale[nrules], norm[nrules]; \
    real gen[NDIM]; \
  } Set

/*********************************************************************/

static void Rule13Alloc(This *t)
{
  static creal w[][nrules] = {
    { .00844923090033615,     .3213775489050763,     .3372900883288987,
     -.8264123822525677,      .6539094339575232 },
    { .023771474018994404,   -.1767341636743844,    -.1644903060344491,
      .306583861409436,      -.2041614154424632},
    { .02940016170142405,     .07347600537466073,    .07707849911634623,
      .002389292538329435,   -.174698151579499 },
    { .006644436465817374,   -.03638022004364754,   -.03804478358506311,
     -.1343024157997222,      .03937939671417803 },
    { .0042536044255016,      .021252979220987123,   .02223559940380806,
      .08833366840533902,     .006974520545933992 },
    { 0,                      .1460984204026913,     .1480693879765931,
      0,                      0 },
    { .0040664827465935255,   .017476132861520992,  4.467143702185815e-6,
      .0009786283074168292,   .0066677021717782585 },
    { .03362231646315497,     .1444954045641582,     .150894476707413,
     -.1319227889147519,      .05512960621544304 },
    { .033200804136503725,    .0001307687976001325, 3.6472001075162155e-5,
      .00799001220015063,     .05443846381278608 },
    { .014093686924979677,    .0005380992313941161,  .000577719899901388,
      .0033917470797606257,   .02310903863953934 },
    { .000977069770327625,    .0001042259576889814,  .0001041757313688177,
      .0022949157182832643,   .01506937747477189 },
    { .007531996943580376,   -.001401152865045733,  -.001452822267047819,
     -.01358584986119197,    -.060570216489018905 },
    { .02577183086722915,     .008041788181514763,   .008338339968783704,
      .04025866859057809,     .04225737654686337},
    { .015625,               -.1420416552759383,    -.147279632923196,
      .003760268580063992,    .02561989142123099 }
  };

  static creal g[] = {
     .12585646717265545,      .3506966822267133,
     .4795480315809981,       .4978005239276064,
     .25,                     .07972723291487795,
     .1904495567970094,       .3291384627633596,
     .43807365825146577,      .499121592026599,
     .4895111329084231,       .32461421628226944,
     .43637106005656195,      .1791307322940614,
     .2833333333333333,       .1038888888888889 };

  enum { nsets = 14, ndim = 2 };

  TYPEDEFSET;

  count n, r;
  Set *first, *last, *s, *x;

  Alloc(first, nsets);
  Clear(first, nsets);

  last = first;
  n = last->n = 1;
  Copy(last->weight, w[0], nrules);

  ++last;
  n += last->n = 2*ndim;
  Copy(last->weight, w[1], nrules);
  last->gen[0] = g[0];

  ++last;
  n += last->n = 2*ndim;
  Copy(last->weight, w[2], nrules);
  last->gen[0] = g[1];

  ++last;
  n += last->n = 2*ndim;
  Copy(last->weight, w[3], nrules);
  last->gen[0] = g[2];

  ++last;
  n += last->n = 2*ndim;
  Copy(last->weight, w[4], nrules);
  last->gen[0] = g[3];

  ++last;
  n += last->n = 2*ndim;
  Copy(last->weight, w[5], nrules);
  last->gen[0] = g[4];

  ++last;
  n += last->n = 2*ndim*(ndim - 1);
  Copy(last->weight, w[6], nrules);
  last->gen[0] = g[5];
  last->gen[1] = g[5];

  ++last;
  n += last->n = 2*ndim*(ndim - 1);
  Copy(last->weight, w[7], nrules);
  last->gen[0] = g[6];
  last->gen[1] = g[6];

  ++last;
  n += last->n = 2*ndim*(ndim - 1);
  Copy(last->weight, w[8], nrules);
  last->gen[0] = g[7];
  last->gen[1] = g[7];

  ++last;
  n += last->n = 2*ndim*(ndim - 1);
  Copy(last->weight, w[9], nrules);
  last->gen[0] = g[8];
  last->gen[1] = g[8];

  ++last;
  n += last->n = 2*ndim*(ndim - 1);
  Copy(last->weight, w[10], nrules);
  last->gen[0] = g[9];
  last->gen[1] = g[9];

  ++last;
  n += last->n = 4*ndim*(ndim - 1);
  Copy(last->weight, w[11], nrules);
  last->gen[0] = g[10];
  last->gen[1] = g[11];

  ++last;
  n += last->n = 4*ndim*(ndim - 1);
  Copy(last->weight, w[12], nrules);
  last->gen[0] = g[12];
  last->gen[1] = g[13];

  ++last;
  n += last->n = 4*ndim*(ndim - 1);
  Copy(last->weight, w[13], nrules);
  last->gen[0] = g[14];
  last->gen[1] = g[15];

  t->rule13.first = first;
  t->rule13.last = last;
  t->rule13.errcoeff[0] = 10;
  t->rule13.errcoeff[1] = 1;
  t->rule13.errcoeff[2] = 5;
  t->rule13.n = n;

  for( s = first; s <= last; ++s )
    for( r = 1; r < nrules - 1; ++r ) {
      creal scale = (s->weight[r] == 0) ? 100 :
        -s->weight[r + 1]/s->weight[r];
      real sum = 0;
      for( x = first; x <= last; ++x )
        sum += x->n*fabs(x->weight[r + 1] + scale*x->weight[r]);
      s->scale[r] = scale;
      s->norm[r] = 1/sum;
    }
}

/*********************************************************************/

static void Rule11Alloc(This *t)
{
  static creal w[][nrules] = {
    { .0009903847688882167,  1.715006248224684,     1.936014978949526,
      .517082819560576,      2.05440450381852 },
    { .0084964717409851,     -.3755893815889209,    -.3673449403754268,
      .01445269144914044,     .013777599884901202 },
    { .00013587331735072814,  .1488632145140549,     .02929778657898176,
     -.3601489663995932,     -.576806291790441 },
    { .022982920777660364,   -.2497046640620823,    -.1151883520260315,
      .3628307003418485,      .03726835047700328 },
    { .004202649722286289,    .1792501419135204,     .05086658220872218,
      .007148802650872729,    .0068148789397772195 },
    { .0012671889041675774,   .0034461267589738897,  .04453911087786469,
     -.09222852896022966,     .057231697338518496 },
    { .0002109560854981544,  -.005140483185555825,  -.022878282571259,
      .01719339732471725,    -.044930187438112855 },
    { .016830857056410086,    .006536017839876424,   .02908926216345833,
     -.102141653746035,       .027292365738663484 },
    { .00021876823557504823, -.00065134549392297,   -.002898884350669207,
     -.007504397861080493,    .000354747395055699 },
    { .009690420479796819,   -.006304672433547204,  -.028059634133074954,
      .01648362537726711,     .01571366799739551 },
    { .030773311284628138,    .01266959399788263,    .05638741361145884,
      .05234610158469334,     .049900992192785674 },
    { .0084974310856038,     -.005454241018647931,  -.02427469611942451,
      .014454323316130661,    .0137791555266677 },
    { .0017749535291258914,   .004826995274768427,   .021483070341828822,
      .003019236275367777,    .0028782064230998723 }
  };

  static creal g[] = {
     .095,                    .25,
     .375,                    .4,
     .4975,                   .49936724991757,
     .38968518428362114,      .49998494965443835,
     .3951318612385894,       .22016983438253684,
     .4774686911397297,       .2189239229503431,
     .4830546566815374,       .2288552938881567 };

  enum { nsets = 13, ndim = 3 };

  TYPEDEFSET;

  count n, r;
  Set *first, *last, *s, *x;

  Alloc(first, nsets);
  Clear(first, nsets);

  last = first;
  n = last->n = 1;
  Copy(last->weight, w[0], nrules);

  ++last;
  n += last->n = 2*ndim;
  Copy(last->weight, w[1], nrules);
  last->gen[0] = g[0];

  ++last;
  n += last->n = 2*ndim;
  Copy(last->weight, w[2], nrules);
  last->gen[0] = g[1];

  ++last;
  n += last->n = 2*ndim;
  Copy(last->weight, w[3], nrules);
  last->gen[0] = g[2];

  ++last;
  n += last->n = 2*ndim;
  Copy(last->weight, w[4], nrules);
  last->gen[0] = g[3];

  ++last;
  n += last->n = 2*ndim;
  Copy(last->weight, w[5], nrules);
  last->gen[0] = g[4];

  ++last;
  n += last->n = 2*ndim*(ndim - 1);
  Copy(last->weight, w[6], nrules);
  last->gen[0] = g[5];
  last->gen[1] = g[5];

  ++last;
  n += last->n = 2*ndim*(ndim - 1);
  Copy(last->weight, w[7], nrules);
  last->gen[0] = g[6];
  last->gen[1] = g[6];

  ++last;
  n += last->n = 4*ndim*(ndim - 1)*(ndim - 2)/3;
  Copy(last->weight, w[8], nrules);
  last->gen[0] = g[7];
  last->gen[1] = g[7];
  last->gen[2] = g[7];

  ++last;
  n += last->n = 4*ndim*(ndim - 1)*(ndim - 2)/3;
  Copy(last->weight, w[9], nrules);
  last->gen[0] = g[8];
  last->gen[1] = g[8];
  last->gen[2] = g[8];

  ++last;
  n += last->n = 4*ndim*(ndim - 1)*(ndim - 2)/3;
  Copy(last->weight, w[10], nrules);
  last->gen[0] = g[9];
  last->gen[1] = g[9];
  last->gen[2] = g[9];

  ++last;
  n += last->n = 4*ndim*(ndim - 1)*(ndim - 2);
  Copy(last->weight, w[11], nrules);
  last->gen[0] = g[10];
  last->gen[1] = g[11];
  last->gen[2] = g[11];

  ++last;
  n += last->n = 4*ndim*(ndim - 1)*(ndim - 2);
  Copy(last->weight, w[12], nrules);
  last->gen[0] = g[12];
  last->gen[1] = g[12];
  last->gen[2] = g[13];

  t->rule11.first = first;
  t->rule11.last = last;
  t->rule11.errcoeff[0] = 4;
  t->rule11.errcoeff[1] = .5;
  t->rule11.errcoeff[2] = 3;
  t->rule11.n = n;

  for( s = first; s <= last; ++s )
    for( r = 1; r < nrules - 1; ++r ) {
      creal scale = (s->weight[r] == 0) ? 100 :
        -s->weight[r + 1]/s->weight[r];
      real sum = 0;
      for( x = first; x <= last; ++x )
        sum += x->n*fabs(x->weight[r + 1] + scale*x->weight[r]);
      s->scale[r] = scale;
      s->norm[r] = 1/sum;
    }
}

/*********************************************************************/

static void Rule9Alloc(This *t)
{
  static creal w[] = {
    -.0023611709677855117884,   .11415390023857325268,
    -.63833920076702389094,     .74849988504685208004,
    -.0014324017033399125142,   .057471507864489725949,
    -.14225104571434243234,    -.062875028738286979989,
     .254591133248959089,     -1.207328566678236261,
     .89567365764160676508,    -.36479356986049146661,
     .0035417564516782676826,  -.072609367395893679605,
     .10557491625218991012,     .0021486025550098687713,
    -.032268563892953949998,    .010636783990231217481,
     .014689102496143490175,    .51134708346467591431,
     .45976448120806344646,     .18239678493024573331,
    -.04508628929435784076,     .21415883524352793401,
    -.027351546526545644722,    .054941067048711234101,
     .11937596202570775297,     .65089519391920250593,
     .14744939829434460168,     .057693384490973483573,
     .034999626602143583822,  -1.3868627719278281436,
    -.2386668732575008879,      .015532417276607053264,
     .0035328099607090870236,   .09231719987444221619,
     .02254314464717892038,     .013675773263272822361,
    -.32544759695960125297,     .0017708782258391338413,
     .0010743012775049343856,   .25150011495314791996 };

  static creal g[] = {
     .47795365790226950619,     .20302858736911986780,
     .44762735462617812882,     .125,
     .34303789878087814570 };

  enum { nsets = 9 };

  TYPEDEFSET;

  ccount ndim = t->ndim;
  ccount twondim = 1 << ndim;
  count dim, n, r;
  Set *first, *last, *s, *x;

  Alloc(first, nsets);
  Clear(first, nsets);

  last = first;
  n = last->n = 1;
  last->weight[0] = ndim*(ndim*(ndim*w[0] + w[1]) + w[2]) + w[3];
  last->weight[1] = ndim*(ndim*(ndim*w[4] + w[5]) + w[6]) - w[7];
  last->weight[2] = ndim*w[8] - last->weight[1];
  last->weight[3] = ndim*(ndim*w[9] + w[10]) - 1 + last->weight[0];
  last->weight[4] = ndim*w[11] + 1 - last->weight[0];

  ++last;
  n += last->n = 2*ndim;
  last->weight[0] = ndim*(ndim*w[12] + w[13]) + w[14];
  last->weight[1] = ndim*(ndim*w[15] + w[16]) + w[17];
  last->weight[2] = w[18] - last->weight[1];
  last->weight[3] = ndim*w[19] + w[20] + last->weight[0];
  last->weight[4] = w[21] - last->weight[0];
  last->gen[0] = g[0];

  ++last;
  n += last->n = 2*ndim;
  last->weight[0] = ndim*w[22] + w[23];
  last->weight[1] = ndim*w[24] + w[25];
  last->weight[2] = w[26] - last->weight[1];
  last->weight[3] = ndim*w[27] + w[28];
  last->weight[4] = -last->weight[0];
  last->gen[0] = g[1];

  ++last;
  n += last->n = 2*ndim;
  last->weight[0] = w[29];
  last->weight[1] = w[30];
  last->weight[2] = -w[29];
  last->weight[3] = w[31];
  last->weight[4] = -w[29];
  last->gen[0] = g[2];

  ++last;
  n += last->n = 2*ndim;
  last->weight[2] = w[32];
  last->gen[0] = g[3];

  ++last;
  n += last->n = 2*ndim*(ndim - 1);
  last->weight[0] = w[33] - ndim*w[12];
  last->weight[1] = w[34] - ndim*w[15];
  last->weight[2] = -last->weight[1];
  last->weight[3] = w[35] + last->weight[0];
  last->weight[4] = -last->weight[0];
  last->gen[0] = g[0];
  last->gen[1] = g[0];

  ++last;
  n += last->n = 4*ndim*(ndim - 1);
  last->weight[0] = w[36];
  last->weight[1] = w[37];
  last->weight[2] = -w[37];
  last->weight[3] = w[38];
  last->weight[4] = -w[36];
  last->gen[0] = g[0];
  last->gen[1] = g[1];

  ++last;
  n += last->n = 4*ndim*(ndim - 1)*(ndim - 2)/3;
  last->weight[0] = w[39];
  last->weight[1] = w[40];
  last->weight[2] = -w[40];
  last->weight[3] = w[39];
  last->weight[4] = -w[39];
  last->gen[0] = g[0];
  last->gen[1] = g[0];
  last->gen[2] = g[0];

  ++last;
  n += last->n = twondim;
  last->weight[0] = w[41]/twondim;
  last->weight[1] = w[7]/twondim;
  last->weight[2] = -last->weight[1];
  last->weight[3] = last->weight[0];
  last->weight[4] = -last->weight[0];
  for( dim = 0; dim < ndim; ++dim )
    last->gen[dim] = g[4];

  t->rule9.first = first;
  t->rule9.last = last;
  t->rule9.errcoeff[0] = 5;
  t->rule9.errcoeff[1] = 1;
  t->rule9.errcoeff[2] = 5;
  t->rule9.n = n;

  for( s = first; s <= last; ++s )
    for( r = 1; r < nrules - 1; ++r ) {
      creal scale = (s->weight[r] == 0) ? 100 :
        -s->weight[r + 1]/s->weight[r];
      real sum = 0;
      for( x = first; x <= last; ++x )
        sum += x->n*fabs(x->weight[r + 1] + scale*x->weight[r]);
      s->scale[r] = scale;
      s->norm[r] = 1/sum;
    }
}

/*********************************************************************/

static void Rule7Alloc(This *t)
{
  static creal w[] = {
     .019417866674748388428,   -.40385257701150182546,
     .64485668767465982223,     .01177982690775806141,
    -.18041318740733609012,    -.088785828081335044443,
     .056328645808285941374,   -.0097089333373741942142,
    -.99129176779582358138,    -.17757165616267008889,
     .12359398032043233572,     .074978148702033690681,
     .55489147051423559776,     .088041241522692771226,
     .021118358455513385083,   -.0099302203239653333087,
    -.064100053285010904179,    .030381729038221007659,
     .0058899134538790307051,  -.0048544666686870971071,
     .35514331232534017777 };

  static creal g[] = {
     .47795365790226950619,     .20302858736911986780,
     .375,                      .34303789878087814570 };

  enum { nsets = 6 };

  TYPEDEFSET;

  ccount ndim = t->ndim;
  ccount twondim = 1 << ndim;
  count dim, n, r;
  Set *first, *last, *s, *x;

  Alloc(first, nsets);
  Clear(first, nsets);

  last = first;
  n = last->n = 1;
  last->weight[0] = ndim*(ndim*w[0] + w[1]) + w[2];
  last->weight[1] = ndim*(ndim*w[3] + w[4]) - w[5];
  last->weight[2] = ndim*w[6] - last->weight[1];
  last->weight[3] = ndim*(ndim*w[7] + w[8]) - w[9];
  last->weight[4] = 1 - last->weight[0];

  ++last;
  n += last->n = 2*ndim;
  last->weight[0] = w[10];
  last->weight[1] = w[11];
  last->weight[2] = -w[10];
  last->weight[3] = w[12];
  last->weight[4] = -w[10];
  last->gen[0] = g[1];

  ++last;
  n += last->n = 2*ndim;
  last->weight[0] = w[13] - ndim*w[0];
  last->weight[1] = w[14] - ndim*w[3];
  last->weight[2] = w[15] - last->weight[1];
  last->weight[3] = w[16] - ndim*w[7];
  last->weight[4] = -last->weight[0];
  last->gen[0] = g[0];

  ++last;
  n += last->n = 2*ndim;
  last->weight[2] = w[17];
  last->gen[0] = g[2];

  ++last;
  n += last->n = 2*ndim*(ndim - 1);
  last->weight[0] = -w[7];
  last->weight[1] = w[18];
  last->weight[2] = -w[18];
  last->weight[3] = w[19];
  last->weight[4] = w[7];
  last->gen[0] = g[0];
  last->gen[1] = g[0];

  ++last;
  n += last->n = twondim;
  last->weight[0] = w[20]/twondim;
  last->weight[1] = w[5]/twondim;
  last->weight[2] = -last->weight[1];
  last->weight[3] = w[9]/twondim;
  last->weight[4] = -last->weight[0];
  for( dim = 0; dim < ndim; ++dim )
    last->gen[dim] = g[3];

  t->rule7.first = first;
  t->rule7.last = last;
  t->rule7.errcoeff[0] = 5;
  t->rule7.errcoeff[1] = 1;
  t->rule7.errcoeff[2] = 5;
  t->rule7.n = n;

  for( s = first; s <= last; ++s )
    for( r = 1; r < nrules - 1; ++r ) {
      creal scale = (s->weight[r] == 0) ? 100 :
        -s->weight[r + 1]/s->weight[r];
      real sum = 0;
      for( x = first; x <= last; ++x )
        sum += x->n*fabs(x->weight[r + 1] + scale*x->weight[r]);
      s->scale[r] = scale;
      s->norm[r] = 1/sum;
    }
}

/*********************************************************************/

static inline void RuleIni(Rule *rule)
{
  rule->first = NULL;
}

/*********************************************************************/

static inline bool RuleIniQ(Rule *rule)
{
  return rule->first == NULL;
}

/*********************************************************************/

static inline void RuleFree(Rule *rule)
{
  free(rule->first);
}

/*********************************************************************/

static real *ExpandFS(cThis *t, cBounds *b, real *g, real *x)
{
  count dim, ndim = t->ndim;

next:
  /* Compute centrally symmetric sum for permutation of G */

  for( dim = 0; dim < ndim; ++dim )
    *x++ = (.5 + g[dim])*b[dim].lower + (.5 - g[dim])*b[dim].upper;

  for( dim = 0; dim < ndim; ) {
    g[dim] = -g[dim];
    if( g[dim++] < 0 ) goto next;
  }

  /* Find next distinct permutation of G and loop back for next sum.
     Permutations are generated in reverse lexicographic order. */

  for( dim = 1; dim < ndim; ++dim ) {
    creal gd = g[dim];
    if( g[dim - 1] > gd ) {
      count i, j = dim, ix = dim, dx = dim - 1;
      for( i = 0; i < --j; ++i ) {
        creal tmp = g[i];
        g[i] = g[j];
        g[j] = tmp;
        if( tmp <= gd ) --dx;
        if( g[i] > gd ) ix = i;
      }
      if( g[dx] <= gd ) dx = ix;
      g[dim] = g[dx];
      g[dx] = gd;
      goto next;
    }
  }

  /* Restore original order to generators */

  for( dim = 0; dim < --ndim; ++dim ) {
    creal tmp = g[dim];
    g[dim] = g[ndim];
    g[ndim] = tmp;
  }

  return x;
}

/*********************************************************************/

static void SampleRule(This *t, ccount iregion)
{
  SAMPLERDEFS;
  TYPEDEFSET;
  Set *first = (Set *)samples->rule->first;
  Set *last = (Set *)samples->rule->last;
  Set *s;
  creal *errcoeff = samples->rule->errcoeff;
  count comp, rul, sn;

  for( s = first; s <= last; ++s )
    if( s->n ) x = ExpandFS(t, b, s->gen, x);

  DoSample(t, n, samples->x, f);

  for( comp = 0; comp < t->ncomp; ++comp ) {
    real sum[nrules];
    creal *f1 = f++;

    Zap(sum);
    for( s = first; s <= last; ++s )
      for( sn = s->n; sn; --sn ) {
        creal fun = *f1;
        f1 += t->ncomp;
        for( rul = 0; rul < nrules; ++rul )
          sum[rul] += fun*s->weight[rul];
      }

    /* Search for the null rule, in the linear space spanned by two
       successive null rules in our sequence, which gives the greatest
       error estimate among all normalized (1-norm) null rules in this
       space. */

    for( rul = 1; rul < nrules - 1; ++rul ) {
      real maxerr = 0;
      for( s = first; s <= last; ++s )
        maxerr = Max(maxerr,
          fabs(sum[rul + 1] + s->scale[rul]*sum[rul])*s->norm[rul]);
      sum[rul] = maxerr;
    }

    r[comp].avg = region->vol*sum[0];
    r[comp].err = region->vol*(
      (errcoeff[0]*sum[1] <= sum[2] && errcoeff[0]*sum[2] <= sum[3]) ?
        errcoeff[1]*sum[1] :
        errcoeff[2]*Max(Max(sum[1], sum[2]), sum[3]) );
  }
}

