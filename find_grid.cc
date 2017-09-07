#include <stdio.h>
#include <vector>
#include <boost/polygon/voronoi.hpp>
#include <assert.h>


// From the boost voronoi-diagram tutorial. Don't ask.
using boost::polygon::voronoi_diagram;
struct Point
{
    int x,y;
    Point(int _x, int _y) : x(_x), y(_y) {}
    Point() {}
};
namespace boost { namespace polygon {
        template <>
        struct geometry_concept<Point> {
            typedef point_concept type;
        };

        template <>
        struct point_traits<Point> {
            typedef int coordinate_type;

            static inline coordinate_type get(const Point& point, orientation_2d orient) {
                return (orient == HORIZONTAL) ? point.x : point.y;
            }
        };
    }}



typedef voronoi_diagram<double> VORONOI;
struct PointDouble
{
    double x,y;
    PointDouble(double _x, double _y) : x(_x), y(_y) {}
    PointDouble() {}
};


#define DEBUG 1



#define Nwant     10
#define SCALE     1000 /* Voronoi diagram is integer-only, so I scale-up to get
                          more resolution */

// tight bound on angle error, loose bound on length error. This is because
// perspective distortion can vary the lengths, but NOT the orientations
#define THRESHOLD_ROW_LENGTH (60*SCALE)
#define THRESHOLD_ROW_COS    0.996 /* 1 degrees */


#define FOR_ALL_ADJACENT_CELLS(c) do {                                  \
    const Point*        pt = &points[c->source_index()];                \
    const VORONOI::edge_type* const e0 = c->incident_edge();            \
    bool first = true;                                                  \
    for(const VORONOI::edge_type* e = e0;                               \
        e0 != NULL && (e != e0 || first);                               \
        e = e->next(), first=false)                                     \
    {                                                                   \
        const VORONOI::cell_type* c_adjacent = e->twin()->cell();       \
                                                                        \
        const Point* pt_adjacent = &points[c_adjacent->source_index()]; \
                                                                        \
        Point delta( pt_adjacent->x - pt->x,                            \
                     pt_adjacent->y - pt->y );

#define FOR_ALL_ADJACENT_CELLS_END() }} while(0)






struct HypothesisStatistics
{
    Point delta;
};

static void fill_initial_hypothesis_statistics(// out
                                               HypothesisStatistics* stats,

                                               // in
                                               const Point* delta0)
{
    stats->delta = *delta0;
}




#define FOR_MATCHING_ADJACENT_CELLS() do {                              \
    HypothesisStatistics stats;                                         \
    fill_initial_hypothesis_statistics(&stats, delta);                  \
    for(int i=0; i<N_remaining; i++)                                    \
    {                                                                   \
        const VORONOI::cell_type* c_adjacent = get_adjacent_cell_along_line(&stats, c, points);


#define FOR_MATCHING_ADJACENT_CELLS_END() \
        c = c_adjacent; }} while(0)




static const VORONOI::cell_type*
get_adjacent_cell_along_line( // out,in.
                              HypothesisStatistics* stats,

                              // in
                              const VORONOI::cell_type* c,
                              const std::vector<Point>& points)
{
    double delta_want_len       = hypot((double)stats->delta.x, (double)stats->delta.y);
    double delta_want_len_recip = 1.0 / delta_want_len;

    FOR_ALL_ADJACENT_CELLS(c)
    {
        double delta_len = hypot( (double)delta.x, (double)delta.y );
        double len_err = delta_want_len - delta_len;
        if( len_err < -THRESHOLD_ROW_LENGTH || THRESHOLD_ROW_LENGTH < len_err )
            continue;

        double cos_err =
            ((double)stats->delta.x * (double)delta.x +
             (double)stats->delta.y * (double)delta.y) *
            delta_want_len_recip / delta_len;
        if( cos_err < THRESHOLD_ROW_COS )
            continue;

#warning returning only one adjacent cell; there could be multiple
        stats->delta = delta;
        return c_adjacent;
    } FOR_ALL_ADJACENT_CELLS_END();

    return NULL;
}

static bool search_along_line( // out
                               PointDouble* delta_mean,

                               // in
                               const Point* delta,
                               const VORONOI::cell_type* c,
                               int N_remaining,

                               const std::vector<Point>& points)
{
    delta_mean->x = (double)delta->x;
    delta_mean->y = (double)delta->y;


    FOR_MATCHING_ADJACENT_CELLS()
    {
        if( c_adjacent == NULL )
            return false;
        delta_mean->x += (double)stats.delta.x;
        delta_mean->y += (double)stats.delta.y;
    }
    FOR_MATCHING_ADJACENT_CELLS_END();

    delta_mean->x /= (double)(N_remaining+1);
    delta_mean->y /= (double)(N_remaining+1);

    return true;
}

static void print_cell_center( const VORONOI::cell_type* c,
                               const std::vector<Point>& points )
{
    const Point* pt = &points[c->source_index()];
    printf("%f %f\n",
           (double)pt->x / (double)SCALE,
           (double)pt->y / (double)SCALE);
}

static void print_along_line( const Point* delta,
                              const VORONOI::cell_type* c,
                              int N_remaining,

                              const std::vector<Point>& points)
{
    FOR_MATCHING_ADJACENT_CELLS()
    {
        print_cell_center(c_adjacent, points);
    } FOR_MATCHING_ADJACENT_CELLS_END();
}

static void dump_interval( const int i_candidate,
                           const int i_pt,
                           const VORONOI::cell_type* c0,
                           const VORONOI::cell_type* c1,
                           const std::vector<Point>& points )
{
    const Point* pt0 = &points[c0->source_index()];
    const Point* pt1 = &points[c1->source_index()];

    double dx = (double)(pt1->x - pt0->x) / (double)SCALE;
    double dy = (double)(pt1->y - pt0->y) / (double)SCALE;
    double length = hypot(dx,dy);
    double angle  = atan2(dy,dx) * 180.0 / M_PI;
    printf("candidate %d point %d, from %f %f to %f %f delta %f %f length %f angle %f\n",
           i_candidate, i_pt,
           (double)pt0->x / (double)SCALE, (double)pt0->y / (double)SCALE,
           (double)pt1->x / (double)SCALE, (double)pt1->y / (double)SCALE,
           dx, dy, length, angle);
}

static void dump_intervals_along_line( int i_candidate,
                                       const Point* delta,
                                       const VORONOI::cell_type* c,
                                       int N_remaining,

                                       const std::vector<Point>& points)
{
    FOR_MATCHING_ADJACENT_CELLS()
    {
        dump_interval(i_candidate, i+1, c, c_adjacent, points);
    } FOR_MATCHING_ADJACENT_CELLS_END();
}


#define CLASSIFICATION_TYPE_LIST(_)             \
    _(UNCLASSIFIED, = 0)                        \
    _(HORIZONTAL,      )                        \
    _(VERTICAL,        )                        \
    _(OUTLIER,         )

#define ENUM_ELEMENT(type,init) type init,
enum ClassificationType
{
 CLASSIFICATION_TYPE_LIST(ENUM_ELEMENT)
};

static const char* type_string(ClassificationType type)
{
#define CASE_RETURN_STRING(type, init) case type: return #type;
    switch(type)
    {
        CLASSIFICATION_TYPE_LIST(CASE_RETURN_STRING)
    default: return "unknown";
    }
}

struct CandidateLine
{
    const VORONOI::cell_type* c0;
    const VORONOI::cell_type* c1;

    PointDouble delta_mean;
    double      spacing_angle;
    double      spacing_length;

    union
    {
        ClassificationType type;
        int bin_index_neg; // if <0, then this is valid
    };
};

double get_spacing_angle( double y, double x )
{
    double angle = 180.0/M_PI * atan2(y,x);
    if( angle < 0 ) angle += 180.0;
    return angle;
}

typedef std::vector<CandidateLine> v_CL;

static bool read_points( std::vector<Point>* points, const char* file )
{
    FILE* fp = fopen(file, "r");
    if( fp == NULL )
    {
        fprintf(stderr, "couldn't open '%s'\n", file);
        return false;
    }

    while(1)
    {
        double x,y;
        int Nread = fscanf(fp, "%lf %lf", &x, &y);
        if(Nread != 2)
            break;

        Point pt( (int)( x * SCALE + 0.5 ),
                  (int)( y * SCALE + 0.5 ) );
        points->push_back(pt);
    }
    fclose(fp);
    return true;
}

static void get_line_candidates( // out
                                 v_CL* line_candidates,

                                 // in
                                 const VORONOI* voronoi,
                                 const std::vector<Point>& points)
{
    for (auto it = voronoi->cells().begin(); it != voronoi->cells().end(); it++ )
    {
        const VORONOI::cell_type* c  = &(*it);

        FOR_ALL_ADJACENT_CELLS(c)
        {
            PointDouble delta_mean;
            if( search_along_line( &delta_mean,
                                   &delta, c_adjacent, Nwant-2, points) )
            {
                double spacing_angle  = get_spacing_angle(delta_mean.y, delta_mean.x);
                double spacing_length = hypot(delta_mean.x, delta_mean.y);

                line_candidates->push_back( CandidateLine({c, c_adjacent, delta_mean,
                                                           spacing_angle, spacing_length}) );

#if defined DEBUG && DEBUG
                printf("pre-classification from %f %f delta_mean %f %f len %f angle %f\n",
                       (double)(pt->x) / (double)SCALE, (double)(pt->y) / (double)SCALE,
                       delta_mean.x / (double)SCALE, delta_mean.y / (double)SCALE,
                       spacing_length / (double)SCALE,
                       spacing_angle);
#endif
            }
        } FOR_ALL_ADJACENT_CELLS_END();
    }
}

struct ClassificationBin
{
    PointDouble delta_mean_sum;
    int N;
};

#define THRESHOLD_BINFIT_LENGTH (120.0*(double)SCALE)
#define THRESHOLD_BINFIT_ANGLE  40.0

static bool fits_in_bin( const CandidateLine* cl,
                         const ClassificationBin* bin,

                         double threshold_binfit_length,
                         double threshold_binfit_angle)
{
    if( bin->N == 0 ) return true;

    double dx         = bin->delta_mean_sum.x/(double)bin->N;
    double dy         = bin->delta_mean_sum.y/(double)bin->N;
    double bin_length = hypot(dx,dy);
    double bin_angle  = get_spacing_angle(dy, dx);

    double length_err = cl->spacing_length - bin_length;
    if(length_err*length_err > threshold_binfit_length*threshold_binfit_length) return false;

    double angle_err = cl->spacing_angle - bin_angle;

    // I want the angular error to be in [-90:90]
    angle_err += 10.0 * 180.0;
    angle_err = remainder(angle_err, 180.0);

    if(angle_err*angle_err > threshold_binfit_angle*threshold_binfit_angle) return false;

    return true;
}

static void push_to_bin( CandidateLine* cl,
                         ClassificationBin* bin,
                         int bin_index )
{
    // I want to accumulate the vector with a consistent absolute direction
    // (angle modu0 180)
    if(bin->delta_mean_sum.x * cl->delta_mean.x +
       bin->delta_mean_sum.y * cl->delta_mean.y >= 0.0 )
    {
        bin->delta_mean_sum.x += cl->delta_mean.x;
        bin->delta_mean_sum.y += cl->delta_mean.y;
    }
    else
    {
        bin->delta_mean_sum.x -= cl->delta_mean.x;
        bin->delta_mean_sum.y -= cl->delta_mean.y;
    }

    cl->bin_index_neg = -bin_index - 1;
    bin->N++;
}

static int gather_unclassified(ClassificationBin* bin,
                               v_CL* line_candidates,
                               int bin_index)
{
    int Nremaining = 0;

    *bin = ClassificationBin({});

    for( auto it = line_candidates->begin(); it != line_candidates->end(); it++ )
    {
        CandidateLine* cl = &(*it);

        if( cl->type != UNCLASSIFIED )
            continue;

        if( fits_in_bin(cl, bin, THRESHOLD_BINFIT_LENGTH, THRESHOLD_BINFIT_ANGLE) )
            push_to_bin(cl, bin, bin_index);
        else
            Nremaining++;
    }

    return Nremaining;
}

static void mark_outliers( v_CL* line_candidates,
                           int bin_index )
{
    for( auto it = line_candidates->begin(); it != line_candidates->end(); it++ )
    {
        CandidateLine* cl = &(*it);

        if( // if we're calling everything remaining and outlier
            bin_index < 0 && cl->type == UNCLASSIFIED ||

            // or we're looking at THIS bin
            cl->bin_index_neg == -bin_index - 1 )
        {
            cl->type = OUTLIER;
        }
    }
}

static void mark_orientation( v_CL* line_candidates,
                              const enum ClassificationType* types )
{
    for( auto it = line_candidates->begin(); it != line_candidates->end(); it++ )
    {
        CandidateLine* cl = &(*it);

        if(      cl->bin_index_neg == -1 ) cl->type = types[0];
        else if( cl->bin_index_neg == -2 ) cl->type = types[1];
    }
}

static bool cluster_line_candidates( v_CL* line_candidates )
{
    // I looked through all my points, and I have candidate sets of points that
    // are
    //
    // - linear-ish
    // - have constant-ish spacing
    //
    // I'm looking for a grid of points. So the spacing angles, lengths should
    // cluster nicely

    // I should have exactly two clusters in spacing_angle/spacing_length space
    // I make an extra bin for outliers
    ClassificationBin bins[3] = {};

    int bin_index = 0;
    while(true)
    {
        ClassificationBin* bin = &bins[bin_index];
        int Nremaining = gather_unclassified( bin, line_candidates, bin_index );

        if( bin->N < Nwant*2 ) // should have enough for both directions
        {
            // this is a bin full of outliers
            mark_outliers( line_candidates, bin_index );
            if( Nremaining == 0)
                // we threw away all the data, and nothing was good.
                return false;
            continue;
        }

        // This was supposedly a good "bin". The last bin is for the outliers,
        // so if it was deemed good, something is off
        if( bin_index >= 2 )
        {
            // too many "good" bins
            return false;
        }

        bin_index++;
        if( Nremaining < Nwant*2 ) // should have enough for both directions
        {
            // only stragglers left. Mark them as outliers and call it good.
            mark_outliers(line_candidates, -1);
            break;
        }
    }

    // alrighty. I now have exactly two bins of exactly the right size. One is
    // "HORIZONTAL" and the other is "VERTICAL". Mark them as such

    // the two orientations should be markedly different from one another. The
    // clustering should have already done this, but I make sure
    enum ClassificationType bin_orientation[2];
    for(int i=0; i<2; i++)
    {
        double angle = get_spacing_angle(bins[i].delta_mean_sum.y, bins[i].delta_mean_sum.x);

        #warning "horizontal/vertical could be ambiguous"
        if( angle > 90-45 && angle < 90+45 )
            bin_orientation[i] = VERTICAL;
        else
            bin_orientation[i] = HORIZONTAL;
    }
    if(bin_orientation[0] == VERTICAL && bin_orientation[1] == VERTICAL)
        return false;
    if(bin_orientation[0] == HORIZONTAL && bin_orientation[1] == HORIZONTAL)
        return false;

    mark_orientation( line_candidates, bin_orientation );
    return true;
}

// Looks through our classification and determines whether things look valid or
// not. Makes no changes to anything
static bool validate_clasification(const v_CL* line_candidates)
{
    // I should have exactly Nwant horizontal lines and Nwant vertical lines
    int Nhorizontal = 0;
    int Nvertical   = 0;

    for( auto it = line_candidates->begin(); it != line_candidates->end(); it++ )
    {
        if(      it->type == HORIZONTAL ) Nhorizontal++;
        else if( it->type == VERTICAL   ) Nvertical++;
    }
    if( Nhorizontal != Nwant ) return false;
    if( Nvertical   != Nwant ) return false;


    // OK then. The horizontal lines should each
#warning complete
    return true;

}

static void dump_candidates( const v_CL* line_candidates,
                             const std::vector<Point>& points )
{
    int N = line_candidates->size();
    for( int i=0; i<N; i++ )
    {
        const CandidateLine* cl = &(*line_candidates)[i];

        dump_interval(i, 0, cl->c0, cl->c1, points);

        const Point* pt0 = &points[cl->c0->source_index()];
        const Point* pt1 = &points[cl->c1->source_index()];

        Point delta({ pt1->x - pt0->x,
                      pt1->y - pt0->y});
        dump_intervals_along_line( i, &delta, cl->c1, Nwant-2, points);
    }
}

static void write_output( const v_CL* line_candidates,
                          const std::vector<Point>& points )
{
    for( auto it = line_candidates->begin(); it != line_candidates->end(); it++ )
    {
        if( it->type == HORIZONTAL )
        {
            print_cell_center(it->c0, points);
            print_cell_center(it->c1, points);

            const Point* pt0 = &points[it->c0->source_index()];
            const Point* pt1 = &points[it->c1->source_index()];

            Point delta({ pt1->x - pt0->x,
                          pt1->y - pt0->y});
            print_along_line( &delta, it->c1, Nwant-2, points);
        }
    }
}

int main(int argc, char* argv[])
{
    if( argc != 2 )
    {
        fprintf(stderr, "need arg on cmdline\n");
        return 1;
    }

    std::vector<Point> points;
    if( !read_points(&points, argv[1]) )
        return 1;

    VORONOI voronoi;
    construct_voronoi(points.begin(), points.end(), &voronoi);

    v_CL line_candidates;
    get_line_candidates(&line_candidates, &voronoi, points);

    if( cluster_line_candidates(&line_candidates))
    {

#if defined DEBUG && DEBUG
        // plot with
        // awk '/post/ {print $3,$13,$4,$6,$7}' | feedgnuplot --dataid --domain --autolegend --square --rangesizeall 3 --with 'vectors size screen 0.01,20 fixed filled'
        for( auto it = line_candidates.begin(); it != line_candidates.end(); it++ )
        {
            const CandidateLine* cl = &(*it);
            const Point*         pt = &points[cl->c0->source_index()];

            printf("post-classification from %f %f delta_mean %f %f len %f angle %f type %s\n",
                   (double)(pt->x) / (double)SCALE, (double)(pt->y) / (double)SCALE,
                   cl->delta_mean.x / (double)SCALE, cl->delta_mean.y / (double)SCALE,
                   cl->spacing_length / (double)SCALE,
                   cl->spacing_angle,
                   type_string(cl->type));
        }
#endif

        if(validate_clasification(&line_candidates))
            write_output(&line_candidates, points);
    }
    return 0;
}
