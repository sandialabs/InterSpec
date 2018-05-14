/*
	partview.cpp
		Partition viewer for Cuba
		last modified 15 Feb 11 th
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <QWidget>
#include <QMainWindow>
#include <QTabWidget>
#include <QAction>
#include <QToolBar>
#include <QPainter>
#include <QPrinter>
#include <QPrintDialog>
#include <QApplication>
#include <QColorGroup>

#include <list>

#include "quit.xpm"
#include "print.xpm"

// --------------------------------------------------------------

class PartitionPlane : public QWidget
{
public:
  PartitionPlane( const int dimx, const int dimy, QWidget *parent = 0 )
  : QWidget(parent), m_dimx(dimx), m_dimy(dimy), m_got(0) {
  }

  void addBound( const int dim, const double lower, const double upper );

  void drawRegion( QPainter *p, const QRect &r );
  void drawRegions( QPainter *p );

  QSize sizeHint() const {
    return QSize(InitialSize, InitialSize);
  }

  QString filename() const {
    return QString("%1-%2.ps").arg(m_dimx).arg(m_dimy);
  }

protected:
  void paintEvent( QPaintEvent * ) {
    QPainter p(this);
    drawRegions( &p );
  }

private:
  int m_dimx, m_dimy, m_got;
  int m_xlower, m_xupper, m_ylower, m_yupper;

  typedef std::list<QRect> regionList;
  typedef regionList::iterator regionIt;
  regionList m_regions;

  enum {
    Hue = 0, // red
    InitialSize = 400,
    CoordScale = 4*InitialSize };
};


void PartitionPlane::addBound( const int dim,
  const double lower, const double upper )
{
  if( dim == m_dimx ) {
    m_xlower = int(CoordScale*lower);
    m_xupper = int(CoordScale*upper);
    m_got |= 1;
  }
  if( dim == m_dimy ) {
    m_ylower = int(CoordScale*lower);
    m_yupper = int(CoordScale*upper);
    m_got |= 2;
  }

  if( m_got == 3 ) {
    m_got = 0;

    const QRect rect = QRect(m_xlower, CoordScale - m_yupper,
      m_xupper - m_xlower + 1, m_yupper - m_ylower + 1);

    const int area = rect.width()*rect.height();
    regionIt r;

    for( r = m_regions.begin(); r != m_regions.end(); ++r ) {
      if( rect == *r ) return;
      if( area > (*r).width()*(*r).height() ) break;
    }
    m_regions.insert(r, rect);

    QPainter p(this);
    drawRegion( &p, rect );
  }
}


void PartitionPlane::drawRegion( QPainter *p, const QRect &r )
{
  p->setWindow(0, 0, CoordScale, CoordScale);

  QColor c;
  const double ratio = r.width()*r.height()/
    double(CoordScale*CoordScale);
  const int saturation = int(255/(M_PI/2)*asin(1 - ratio));
  c.setHsv( Hue, saturation, 255 );
  p->setBrush(c);

  p->setPen(palette().foreground().color());

  p->drawRect(r);
}


void PartitionPlane::drawRegions( QPainter *p )
{
  for( regionIt r = m_regions.begin(); r != m_regions.end(); ++r )
    drawRegion( p, *r );
}


// --------------------------------------------------------------

class PartitionViewer : public QMainWindow
{
  Q_OBJECT

public:
  PartitionViewer( QWidget *parent );
  void addPlane( const int dimx, const int dimy );
  void addBound( const int dim, const double lower, const double upper );
  int count() const { return m_tabs->count(); }
  void tabupdate() { 
    if( m_tabs->currentWidget() != 0 ) m_tabs->currentWidget()->update();
  }

public slots:
  void print();

private:
  QTabWidget *m_tabs;
  QPrinter *m_printer;
};


PartitionViewer::PartitionViewer( QWidget *parent = 0 )
: QMainWindow(parent)
{
  setWindowTitle(tr("Cuba Partition Viewer"));

  QToolBar *toolbar = new QToolBar(this);
  addToolBar(Qt::LeftToolBarArea, toolbar);

  QAction *quit = new QAction( QPixmap(quit_xpm), tr("&Quit"), this );
  connect( quit, SIGNAL(activated()), qApp, SLOT(quit()) );
  toolbar->addAction(quit);

#ifndef QT_NO_PRINTER
  QAction *print = new QAction( QPixmap(print_xpm), tr("&Print..."), this );
  connect( print, SIGNAL(activated()), this, SLOT(print()) );
  toolbar->addAction(print);

  m_printer = new QPrinter;
#endif

  m_tabs = new QTabWidget(this);
  setCentralWidget(m_tabs);
}


void PartitionViewer::addPlane( const int dimx, const int dimy )
{
  PartitionPlane *plane = new PartitionPlane(dimx, dimy, this);
  m_tabs->addTab( plane, tr("%1-%2 plane").arg(dimx).arg(dimy) );
}


void PartitionViewer::addBound( const int dim,
  const double lower, const double upper )
{
  for( int index = 0; index < m_tabs->count(); ++index ) {
    PartitionPlane *plane = (PartitionPlane *)m_tabs->widget(index);
    if( plane ) plane->addBound(dim, lower, upper);
  }
}


void PartitionViewer::print()
{
#ifndef QT_NO_PRINTER
  PartitionPlane *plane = (PartitionPlane *)m_tabs->currentWidget();
  if( !plane ) return;

  QPrintDialog printDialog(m_printer, this);
  if( printDialog.exec() == QDialog::Accepted ) {
    QPainter p;
    if( p.begin(m_printer) ) {
      p.setViewport( QRect(QPoint(0, 0), plane->sizeHint()) );
      plane->drawRegions( &p );
    }
  }
#endif
}


#include "partview.moc"

// --------------------------------------------------------------

int main( int argc, char **argv )
{
  QApplication app(argc, argv);
  PartitionViewer partview;

  argc = (argc - 1) & -2;

  for( int arg = 0; arg < argc; ) {
    const int dimx = atoi(argv[++arg]);
    const int dimy = atoi(argv[++arg]);
    if( dimx > 0 && dimy > 0 ) partview.addPlane(dimx, dimy);
  }

  if( partview.count() == 0 ) {
    fprintf(stderr, "Usage:  %s dimx dimy ...\n"
      "reads Cuba's verbose = 3 output from stdin and displays\n"
      "the dimx-dimy plane of the tessellation on screen.\n"
      "Each pair of dimensions is shown in a separate window.\n\n",
      argv[0]);
    exit(1);
  }

  partview.show();

  int dim = 0;
  char line[128];

  while( fgets(line, sizeof line, stdin) ) {
    double lower, upper;

    fputs(line, stdout);

    if( sscanf(line, "%*[^(](%lf) - (%lf)", &lower, &upper) == 2 )
      partview.addBound(++dim, lower, upper);
    else dim = 0;

    app.processEvents();
  }

  fflush(stdout);
  partview.tabupdate();
  return app.exec();
}

