#include <sys/stat.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <set>
#include <climits>
#include <boost/polygon/voronoi.hpp>
#include <assert.h>
#include "point.hh"
#include "mrgingham.hh"
#include "mrgingham-internal.h"

using namespace mrgingham;


// From the boost voronoi-diagram tutorial. Don't ask.
using boost::polygon::voronoi_diagram;
namespace boost { namespace polygon {
        template <>
        struct geometry_concept<PointInt> {
            typedef point_concept type;
        };

        template <>
        struct point_traits<PointInt> {
            typedef int coordinate_type;

            static inline coordinate_type get(const PointInt& point, orientation_2d orient) {
                return (orient == HORIZONTAL) ? point.x : point.y;
            }
        };
    }}



typedef voronoi_diagram<double> VORONOI;


// hard-coding 10x10 grids
#define Nwant     10


/*
Voronoi library terminology. Each "cell" c has an edge = c->incident_edge(). The
edge e points FROM the cell e->cell(). The edges are a linked list. e->next() is
the next edge from c, moving clockwise. e->prev() moves counterclockwise(). An
edge e A->B has an inverse e->twin() B->A.

For any voronoi vertex ("cell" in the terminology of this voronoi library) I
want look at all the neighboring cells. The simplest thing to do is to look at
all the vertices directly connected to this vertex with an edge. This works most
of the time. But for very skewed views of a chessboard this sometimes isn't
enough, and I look at the in-between vertices too. From vertex A below I would
previously consider E and B and C as neighbors, but I now also want to consider
D

E--B-----D
| / \   /
|/   \ /
A-----C

How do I find D? Let e be A->B. Then D =
e->twin()->prev()->prev()->twin()->cell()

There're a few topological wrinkles that this scheme creates. First off, at the
boundary of the graph, there may not exist an "inbetween" vertex. In the example
above there's no vertex between C and E from A (looking clockwise). I check by
making sure that e->twin()->prev()->twin()->cell() == e->next()->twin()->cell().
In the "good" example of D being between B and C, both should point to C.

I can also have more complex topologies, which would create duplicate neighbors.
For instance if I have this:

    ---F
  -/  /|
 /   / |
A---H  |
 \   \ |
  -\  \|
    ---G

Looking from A, vertex G lies between F,H. But G is directly connected to A, so
I would already use this as a neighbor. There's a whole class of such issues,
and I avoid them by only selecting vertices that are monotonic in angle, moving
around A. So in the above example the midpoint MUST lie to the right of F and to
the left of H. But G lies to the right of H, so I do not report an in-between
vertex between F and H.

This is all a bit heuristic, but is sufficient for the purposes of chessboard
grid finding.

 */
#define FOR_ALL_ADJACENT_CELLS(c) do {                                  \
    const PointInt*        pt = &points[c->source_index()];             \
    const VORONOI::edge_type* const __e0 = c->incident_edge();          \
    bool __first = true;                                                \
    for(const VORONOI::edge_type* __e = __e0;                           \
        __e0 != NULL && (__e != __e0 || __first);                       \
        __e = __e->next(), __first=false)                               \
    {                                                                   \
        /* Look for two neighbors using this edge. The edge itself, */  \
        /* and the vertex between this edge and the next edge */        \
        for(int __i=0; __i<2; __i++)                                    \
        {                                                               \
            const VORONOI::cell_type* c_adjacent = __e->twin()->cell(); \
            const PointInt* pt_adjacent = &points[c_adjacent->source_index()]; \
            if(__i == 0)                                                \
                ; /* use c_adjacent, pt_adjacent above */               \
            else                                                        \
            {                                                           \
                const VORONOI::cell_type* c0 = c_adjacent;              \
                const PointInt* pt0 = pt_adjacent;                      \
                const VORONOI::cell_type* c1 = __e->next()->twin()->cell(); \
                const PointInt* pt1 = &points[c1->source_index()];      \
                PointInt v0( pt0->x - pt->x, pt0->y - pt->y );          \
                PointInt v1( pt1->x - pt->x, pt1->y - pt->y );          \
                                                                        \
                /* check for out-of-bounds midpoints in two ways */     \
                /* just one method is probably enough, but I want to make sure */ \
                                                                        \
                /* if this edge and the next edge don't form an acute angle, */ \
                /* then we're at a graph boundary, and there's no midpoint */ \
                if((long)v1.x*(long)v0.y > (long)v0.x*(long)v1.y) \
                    continue;                                           \
                                                                        \
                /* if this and the next edge don't form a triangle, */ \
                /* then we're at a graph boundary, and there's no midpoint */ \
                if(__e->twin()->prev()->twin()->cell() != __e->next()->twin()->cell() ) \
                    continue;                                           \
                                                                        \
                /* get the midpoint vertex. It must lie angularly between its two neighbors */ \
                c_adjacent = __e->twin()->prev()->prev()->twin()->cell(); \
                const PointInt* ptmid = &points[c_adjacent->source_index()]; \
                PointInt vmid( ptmid->x - pt->x, ptmid->y - pt->y );    \
                if((long)v1  .x*(long)vmid.y > (long)vmid.x*(long)v1  .y) \
                    continue;                                           \
                if((long)vmid.x*(long)v0  .y > (long)v0  .x*(long)vmid.y) \
                    continue;                                           \
                pt_adjacent = &points[c_adjacent->source_index()];      \
            }                                                           \
                                                                        \
            PointInt delta( pt_adjacent->x - pt->x,                     \
                            pt_adjacent->y - pt->y );

#define FOR_ALL_ADJACENT_CELLS_END() }}} while(0)




struct CandidateSequence
{
    // First two cells and the last cell. The rest of the sequence is
    // constructed by following the best path from these two cells. The rules
    // that define this "best" path are consistent, so we don't store the path
    // itself, but recompute it each time it is needed
    const VORONOI::cell_type* c0;
    const VORONOI::cell_type* c1;
    const VORONOI::cell_type* clast;

    PointDouble delta_mean;
    double      spacing_angle;
    double      spacing_length;
};




struct HypothesisStatistics
{
    PointInt  delta_last;

    double length_ratio_sum;
    int    length_ratio_N;
};

static void fill_initial_hypothesis_statistics(// out
                                               HypothesisStatistics* stats,

                                               // in
                                               const PointInt* delta0)
{
    stats->delta_last       = *delta0;
    stats->length_ratio_sum = 0.0;
    stats->length_ratio_N   = 0;
}



// need delta, N_remaining, c, points
#define FOR_MATCHING_ADJACENT_CELLS(debug_sequence_pointscale) do {     \
    HypothesisStatistics stats;                                         \
    fill_initial_hypothesis_statistics(&stats, delta);                  \
    for(int i=0; i<N_remaining; i++)                                    \
    {                                                                   \
        const VORONOI::cell_type* c_adjacent = get_adjacent_cell_along_sequence(&stats, c, points, debug_sequence_pointscale);


#define FOR_MATCHING_ADJACENT_CELLS_END() \
        c = c_adjacent; }} while(0)




// tight bound on angle error, loose bound on length error. This is because
// perspective distortion can vary the lengths, but NOT the orientations
#define THRESHOLD_SPACING_LENGTH                 (80.*FIND_GRID_SCALE)
#define THRESHOLD_SPACING_COS                    0.984 /* 10 degrees */
#define THRESHOLD_SPACING_LENGTH_RATIO_MIN       0.7
#define THRESHOLD_SPACING_LENGTH_RATIO_MAX       1.4
#define THRESHOLD_SPACING_LENGTH_RATIO_DEVIATION 0.25

static const VORONOI::cell_type*
get_adjacent_cell_along_sequence( // out,in.
                                 HypothesisStatistics* stats,

                                 // in
                                 const VORONOI::cell_type* c,
                                 const std::vector<PointInt>& points,
                                 int debug_sequence_pointscale /* <=0 means "no debugging" */ )
{
    // We're given a voronoi cell, and some properties that a potential next
    // cell in the sequence should match. I look through all the voronoi
    // neighbors of THIS cell, and return the first one that matches all my
    // requirements. Multiple neighboring cells COULD match, but I'm assuming
    // clean data, so this possibility is ignored.
    //
    // A matching next cell should have the following properties, all
    // established by the previous cells in the sequence. The next cell should
    // match all of these:
    //
    // - located along an expected direction (tight bound on angle)
    //
    // - should be an expected distance away (loose bound on absolute distance)
    //
    // - the ratio of successive distances in a sequence should be constant-ish
    //
    // The points in a sequence come from projections of equidistantly-spaced
    // points in a straight line, so the projected DIRECTIONS should match very
    // well. The distances may not, if the observed object is tilted and is
    // relatively close to the camera, but each successive distance will vary ~
    // geometrically due to perspective effects, or it the distances will all be
    // roughly constant, which is still geometric, technically

    const PointInt& delta_last = stats->delta_last;

    double delta_last_length = hypot((double)delta_last.x, (double)delta_last.y);

    FOR_ALL_ADJACENT_CELLS(c)
    {
        if(debug_sequence_pointscale > 0)
            fprintf(stderr, "Considering connection in sequence from (%d,%d) -> (%d,%d); delta (%d,%d) ..... \n",
                    pt->x          / debug_sequence_pointscale,
                    pt->y          / debug_sequence_pointscale,
                    pt_adjacent->x / debug_sequence_pointscale,
                    pt_adjacent->y / debug_sequence_pointscale,
                    delta.x        / debug_sequence_pointscale,
                    delta.y        / debug_sequence_pointscale);

        double delta_length = hypot( (double)delta.x, (double)delta.y );

        double cos_err =
            ((double)delta_last.x * (double)delta.x +
             (double)delta_last.y * (double)delta.y) /
            (delta_last_length * delta_length);
        if( cos_err < THRESHOLD_SPACING_COS )
        {
            if(debug_sequence_pointscale > 0)
                fprintf(stderr, "..... rejecting. Angle is wrong. I wanted cos_err>=threshold, but saw %f<%f\n",
                        cos_err, THRESHOLD_SPACING_COS);
            continue;
        }

        double length_err = delta_last_length - delta_length;
        if( length_err < -THRESHOLD_SPACING_LENGTH ||
            length_err >  THRESHOLD_SPACING_LENGTH )
        {
            if(debug_sequence_pointscale > 0)
                fprintf(stderr, "..... rejecting. Lengths are wrong. I wanted abs(length_err)<=threshold, but saw %f>%f\n",
                        fabs(length_err), THRESHOLD_SPACING_LENGTH);
            continue;
        }

        double length_ratio = delta_length / delta_last_length;
        if( length_ratio < THRESHOLD_SPACING_LENGTH_RATIO_MIN ||
            length_ratio > THRESHOLD_SPACING_LENGTH_RATIO_MAX )
        {
            if(debug_sequence_pointscale > 0)
                fprintf(stderr, "..... rejecting. Lengths are wrong. I wanted abs(length_ratio)<=threshold, but saw %f<%f or %f>%f\n",
                        length_ratio, THRESHOLD_SPACING_LENGTH_RATIO_MIN,
                        length_ratio, THRESHOLD_SPACING_LENGTH_RATIO_MAX);
            continue;
        }

        // I compute the mean and look at the deviation from the CURRENT mean. I
        // ignore the first few points, since the mean is unstable then. This is
        // OK, however, since I'm going to find and analyze the same sequence in
        // the reverse order, and this will cover the other end
        if( stats->length_ratio_N > 2 )
        {
            double length_ratio_mean = stats->length_ratio_sum / (double)stats->length_ratio_N;

            double length_ratio_deviation = length_ratio - length_ratio_mean;
            if( length_ratio_deviation < -THRESHOLD_SPACING_LENGTH_RATIO_DEVIATION ||
                length_ratio_deviation >  THRESHOLD_SPACING_LENGTH_RATIO_DEVIATION )
            {
                if(debug_sequence_pointscale > 0)
                    fprintf(stderr, "..... rejecting. Lengths are wrong. I wanted abs(length_ratio_deviation)<=threshold, but saw %f>%f\n",
                            fabs(length_ratio_deviation), THRESHOLD_SPACING_LENGTH_RATIO_DEVIATION);
                continue;
            }
        }

        stats->length_ratio_sum += length_ratio;
        stats->length_ratio_N++;

        stats->delta_last        = delta;

        if(debug_sequence_pointscale > 0)
            fprintf(stderr, "..... accepting!\n\n");
        return c_adjacent;

    } FOR_ALL_ADJACENT_CELLS_END();

    return NULL;
}

static
const VORONOI::cell_type* search_along_sequence( // out
                                                PointDouble* delta_mean,

                                                // in
                                                const PointInt* delta,
                                                const VORONOI::cell_type* c,
                                                int N_remaining,

                                                const std::vector<PointInt>& points,
                                                int debug_sequence_pointscale )
{
    delta_mean->x = (double)delta->x;
    delta_mean->y = (double)delta->y;

    const VORONOI::cell_type* clast = NULL;
    FOR_MATCHING_ADJACENT_CELLS(debug_sequence_pointscale)
    {
        if( c_adjacent == NULL )
            return NULL;
        delta_mean->x += (double)stats.delta_last.x;
        delta_mean->y += (double)stats.delta_last.y;

        if(i == N_remaining-1)
            clast = c_adjacent;
    }
    FOR_MATCHING_ADJACENT_CELLS_END();

    delta_mean->x /= (double)(N_remaining+1);
    delta_mean->y /= (double)(N_remaining+1);

    return clast;
}

static void output_point( std::vector<PointDouble>& points_out,
                          const VORONOI::cell_type* c,
                          const std::vector<PointInt>& points )
{
    const PointInt* pt = &points[c->source_index()];
    points_out.push_back( PointDouble( (double)pt->x / (double)FIND_GRID_SCALE,
                                       (double)pt->y / (double)FIND_GRID_SCALE) );
}

static void output_points_along_sequence( std::vector<PointDouble>& points_out,
                                          const PointInt* delta,
                                          const VORONOI::cell_type* c,
                                          int N_remaining,

                                          const std::vector<PointInt>& points)
{
    FOR_MATCHING_ADJACENT_CELLS(-1)
    {
        output_point(points_out, c_adjacent, points);
    } FOR_MATCHING_ADJACENT_CELLS_END();
}

static void output_row( std::vector<PointDouble>& points_out,
                        const CandidateSequence& row,
                        const std::vector<PointInt>& points )
{
    output_point(points_out, row.c0, points);
    output_point(points_out, row.c1, points);

    const PointInt* pt0 = &points[row.c0->source_index()];
    const PointInt* pt1 = &points[row.c1->source_index()];

    PointInt delta({ pt1->x - pt0->x,
                     pt1->y - pt0->y});
    output_points_along_sequence( points_out, &delta, row.c1, Nwant-2, points);
}

// dumps the voronoi diagram to a self-plotting vnlog
#define DUMP_FILENAME_VORONOI "/tmp/mrgingham-2-voronoi.vnl"
static void dump_voronoi( const VORONOI* voronoi,
                          const std::vector<PointInt>& points )
{
    FILE* fp = fopen(DUMP_FILENAME_VORONOI, "w");
    assert(fp);

    // the kernel limits the #! line to 127 characters, so I abbreviate
    fprintf(fp, "#!/usr/bin/feedgnuplot --domain --dataid --with 'lines linecolor 0' --square --maxcurves 100000 --set 'yrange [:] rev'\n");
    fprintf(fp, "# x id_edge y\n");

    int i_edge = 0;
    for (auto it = voronoi->cells().begin(); it != voronoi->cells().end(); it++ )
    {
        const VORONOI::cell_type* c = &(*it);

        const PointInt*        pt0 = &points[c->source_index()];
        const VORONOI::edge_type* const e0 = c->incident_edge();
        bool first = true;
        for(const VORONOI::edge_type* e = e0;
            e0 != NULL && (e != e0 || first);
            e = e->next(), first=false)
        {
            const VORONOI::cell_type* c_adjacent = e->twin()->cell();
            const PointInt* pt1 = &points[c_adjacent->source_index()];
            fprintf(fp, "%f %d %f\n", pt0->x/(double)FIND_GRID_SCALE, i_edge, pt0->y/(double)FIND_GRID_SCALE);
            fprintf(fp, "%f %d %f\n", pt1->x/(double)FIND_GRID_SCALE, i_edge, pt1->y/(double)FIND_GRID_SCALE);
            i_edge++;
        }
    }
    fclose(fp);
    chmod(DUMP_FILENAME_VORONOI,
          S_IRUSR | S_IRGRP | S_IROTH |
          S_IWUSR | S_IWGRP |
          S_IXUSR | S_IXGRP | S_IXOTH);
    fprintf(stderr, "Wrote self-plotting voronoi diagram to " DUMP_FILENAME_VORONOI "\n");
}

static void dump_interval( FILE* fp,
                           const int i_candidate,
                           const int i_pt,
                           const VORONOI::cell_type* c0,
                           const VORONOI::cell_type* c1,
                           const std::vector<PointInt>& points )
{
    const PointInt* pt0 = &points[c0->source_index()];

    if( c1 == NULL )
    {
        fprintf(fp,
                "%d %d %f %f - - - - - -\n",
                i_candidate, i_pt,
                (double)pt0->x / (double)FIND_GRID_SCALE, (double)pt0->y / (double)FIND_GRID_SCALE);
        return;
    }

    const PointInt* pt1 = &points[c1->source_index()];
    double dx = (double)(pt1->x - pt0->x) / (double)FIND_GRID_SCALE;
    double dy = (double)(pt1->y - pt0->y) / (double)FIND_GRID_SCALE;
    double length = hypot(dx,dy);
    double angle  = atan2(dy,dx) * 180.0 / M_PI;
    fprintf(fp,
            "%d %d %f %f %f %f %f %f %f %f\n",
           i_candidate, i_pt,
           (double)pt0->x / (double)FIND_GRID_SCALE, (double)pt0->y / (double)FIND_GRID_SCALE,
           (double)pt1->x / (double)FIND_GRID_SCALE, (double)pt1->y / (double)FIND_GRID_SCALE,
           dx, dy, length, angle);
}

static void dump_intervals_along_sequence( FILE* fp,
                                           const CandidateSequence* cs,
                                           int i_candidate,

                                           const std::vector<PointInt>& points)
{
    int N_remaining = Nwant-1;

    dump_interval(fp, i_candidate, 0, cs->c0, cs->c1, points);

    const PointInt* pt0 = &points[cs->c0->source_index()];
    const PointInt* pt1 = &points[cs->c1->source_index()];

    PointInt _delta({ pt1->x - pt0->x,
                      pt1->y - pt0->y});
    const PointInt* delta = &_delta;

    const VORONOI::cell_type* c = cs->c1;

    FOR_MATCHING_ADJACENT_CELLS(-1)
    {
        dump_interval(fp, i_candidate, i+1, c,
                      i+1 == Nwant-1 ? NULL : c_adjacent,
                      points);
    } FOR_MATCHING_ADJACENT_CELLS_END();
}

static
double get_spacing_angle( double y, double x )
{
    double angle = 180.0/M_PI * atan2(y,x);
    if( angle < 0 ) angle += 180.0;
    return angle;
}

typedef std::vector<CandidateSequence>  v_CS;

typedef std::map<unsigned int, std::vector<int> > SequenceIndicesFromPoint;

struct outer_cycle
{
    short e[4];
};


static void get_sequence_candidates( // out
                                     v_CS* sequence_candidates,

                                     // in
                                     const VORONOI* voronoi,
                                     const std::vector<PointInt>& points,

                                     // for debugging
                                     const debug_sequence_t& debug_sequence)
{
    const VORONOI::cell_type* tracing_c = NULL;

    int debug_sequence_pointscale = -1;
    if(debug_sequence.dodebug)
    {
        // we're tracing some point. I find the nearest voronoi vertex, and
        // debug_sequence that
        unsigned long d2 = (unsigned long)(-1L); // max at first
        debug_sequence_pointscale = FIND_GRID_SCALE;
        for (auto it = voronoi->cells().begin(); it != voronoi->cells().end(); it++ )
        {
            const VORONOI::cell_type* c  = &(*it);
            const PointInt*           pt = &points[c->source_index()];
            long dx = (long)(pt->x - debug_sequence_pointscale*debug_sequence.pt.x);
            long dy = (long)(pt->y - debug_sequence_pointscale*debug_sequence.pt.y);
            unsigned long d2_here = (unsigned long)(dx*dx + dy*dy);
            if(d2_here < d2)
            {
                d2 = d2_here;
                tracing_c = c;
            }
        }
        const PointInt* pt = &points[tracing_c->source_index()];
        fprintf(stderr, "============== Looking at sequences from (%d,%d)\n",
                pt->x / debug_sequence_pointscale,
                pt->y / debug_sequence_pointscale);
    }

    for (auto it = voronoi->cells().begin(); it != voronoi->cells().end(); it++ )
    {
        const VORONOI::cell_type* c  = &(*it);

        bool debug_sequence = ( c == tracing_c );

        FOR_ALL_ADJACENT_CELLS(c)
        {
            if(c == tracing_c)
                fprintf(stderr, "\n\n====== Looking at adjacent point (%d,%d)\n",
                        pt_adjacent->x / debug_sequence_pointscale,
                        pt_adjacent->y / debug_sequence_pointscale);

            PointDouble delta_mean;
            const VORONOI::cell_type* clast =
                search_along_sequence( &delta_mean,
                                       &delta, c_adjacent, Nwant-2, points,
                                       (c == tracing_c) ? debug_sequence_pointscale : -1 );
            if( clast )
            {
                double spacing_angle  = get_spacing_angle(delta_mean.y, delta_mean.x);
                double spacing_length = hypot(delta_mean.x, delta_mean.y);

                sequence_candidates->push_back( CandidateSequence({c, c_adjacent, clast, delta_mean,
                                                                   spacing_angle, spacing_length}) );
            }
        } FOR_ALL_ADJACENT_CELLS_END();
    }
}

static void get_candidate_point( unsigned int* cs_point,
                                 const VORONOI::cell_type* c )
{
    *cs_point = c->source_index();
}
static void get_candidate_points_along_sequence( unsigned int* cs_points,

                                                 const PointInt* delta,
                                                 const VORONOI::cell_type* c,
                                                 int N_remaining,

                                                 const std::vector<PointInt>& points)
{
    FOR_MATCHING_ADJACENT_CELLS(-1)
    {
        get_candidate_point(cs_points, c_adjacent);
        cs_points++;
    } FOR_MATCHING_ADJACENT_CELLS_END();
}
static void get_candidate_points( unsigned int* cs_points,
                                  const CandidateSequence* cs,
                                  const std::vector<PointInt>& points )
{
    get_candidate_point( &cs_points[0], cs->c0 );
    get_candidate_point( &cs_points[1], cs->c1 );

    const PointInt* pt0 = &points[cs->c0->source_index()];
    const PointInt* pt1 = &points[cs->c1->source_index()];

    PointInt delta({ pt1->x - pt0->x,
                  pt1->y - pt0->y});
    get_candidate_points_along_sequence(&cs_points[2], &delta, cs->c1, Nwant-2, points);
}


// dumps a terse self-plotting vnlog visualization of sequence candidates, and a
// more detailed vnlog containing more data
#define DUMP_BASENAME_ALL_SEQUENCE_CANDIDATES     "/tmp/mrgingham-3-candidates"
#define DUMP_BASENAME_OUTER_EDGES                 "/tmp/mrgingham-4-outer-edges"
#define DUMP_BASENAME_OUTER_EDGE_CYCLES           "/tmp/mrgingham-5-outer-edge-cycles"
#define DUMP_BASENAME_IDENTIFIED_OUTER_EDGE_CYCLE "/tmp/mrgingham-6-identified-outer-edge-cycle"
#define dump_candidates(basename, sequence_candidates, outer_edges, points) \
    _dump_candidates( basename".vnl", basename"-detailed.vnl",  \
                      sequence_candidates, outer_edges, points )
static void _dump_candidates(const char* filename_sparse,
                             const char* filename_dense,
                             const v_CS* sequence_candidates,
                             const std::vector<int>* outer_edges,
                             const std::vector<PointInt>& points)
{
    FILE* fp = fopen(filename_sparse, "w");
    assert(fp);

    // the kernel limits the #! line to 127 characters, so I abbreviate
    fprintf(fp, "#!/usr/bin/feedgnuplot --dom --aut --square --rangesizea 3 --w 'vec size screen 0.01,20 fixed fill' --set 'yr [:] rev'\n");
    fprintf(fp, "# fromx type fromy deltax deltay\n");

    if(outer_edges == NULL)
        for( auto it = sequence_candidates->begin(); it != sequence_candidates->end(); it++ )
        {
            const CandidateSequence* cs = &(*it);
            const PointInt*          pt = &points[cs->c0->source_index()];

            fprintf(fp,
                    "%f %f %f %f\n",
                    (double)(pt->x) / (double)FIND_GRID_SCALE,
                    (double)(pt->y) / (double)FIND_GRID_SCALE,
                    cs->delta_mean.x / (double)FIND_GRID_SCALE,
                    cs->delta_mean.y / (double)FIND_GRID_SCALE);
        }
    else
        for( auto it = outer_edges->begin(); it != outer_edges->end(); it++ )
        {
            const CandidateSequence* cs = &((*sequence_candidates)[*it]);
            const PointInt*          pt = &points[cs->c0->source_index()];

            fprintf(fp,
                    "%f %f %f %f\n",
                    (double)(pt->x) / (double)FIND_GRID_SCALE,
                    (double)(pt->y) / (double)FIND_GRID_SCALE,
                    cs->delta_mean.x / (double)FIND_GRID_SCALE,
                    cs->delta_mean.y / (double)FIND_GRID_SCALE);
        }

    fclose(fp);
    chmod(filename_sparse,
          S_IRUSR | S_IRGRP | S_IROTH |
          S_IWUSR | S_IWGRP |
          S_IXUSR | S_IXGRP | S_IXOTH);
    fprintf(stderr, "Wrote self-plotting sequence-candidate dump to %s\n",
            filename_sparse);


    // detailed
    fp = fopen(filename_dense, "w");
    assert(fp);

    fprintf(fp, "# candidateid pointid fromx fromy tox toy deltax deltay len angle\n");

    if(outer_edges == NULL)
    {
        int N = sequence_candidates->size();
        for( int i=0; i<N; i++ )
            dump_intervals_along_sequence( fp,
                                           &(*sequence_candidates)[i],
                                           i, points);
    }
    else
    {
        int N = outer_edges->size();
        for( int i=0; i<N; i++ )
            dump_intervals_along_sequence( fp,
                                           &(*sequence_candidates)[(*outer_edges)[i]],
                                           i, points);
    }
    fclose(fp);
    fprintf(stderr, "Wrote detailed sequence-candidate dump to %s\n",
            filename_dense);
}

static void dump_outer_edge_cycles(const std::vector<outer_cycle>& outer_cycles,
                                   const std::vector<int>&         outer_edges,
                                   const v_CS&                     sequence_candidates,
                                   const std::vector<PointInt>&    points)
{
    FILE* fp = fopen(DUMP_BASENAME_OUTER_EDGE_CYCLES, "w");
    assert(fp);

    // the kernel limits the #! line to 127 characters, so I abbreviate
    fprintf(fp, "#!/usr/bin/feedgnuplot --datai --dom --aut --square --rangesizea 3 --w 'vec size screen 0.01,20 fixed fill' --set 'yr [:] rev'\n");
    fprintf(fp, "# fromx type fromy deltax deltay\n");

    for( int i_cycle=0; i_cycle<(int)outer_cycles.size(); i_cycle++ )
    {
        for(int i_edge = 0; i_edge<4; i_edge++ )
        {
            const CandidateSequence* cs = &sequence_candidates[outer_edges[ outer_cycles[i_cycle].e[i_edge] ]];
            const PointInt*          pt = &points[cs->c0->source_index()];

            fprintf(fp,
                    "%f %d %f %f %f\n",
                    (double)(pt->x) / (double)FIND_GRID_SCALE,
                    i_cycle,
                    (double)(pt->y) / (double)FIND_GRID_SCALE,
                    cs->delta_mean.x / (double)FIND_GRID_SCALE,
                    cs->delta_mean.y / (double)FIND_GRID_SCALE);
        }
    }

    fclose(fp);
    chmod(DUMP_BASENAME_OUTER_EDGE_CYCLES,
          S_IRUSR | S_IRGRP | S_IROTH |
          S_IWUSR | S_IWGRP |
          S_IXUSR | S_IXGRP | S_IXOTH);
    fprintf(stderr, "Wrote outer edge cycle dump to %s\n",
            DUMP_BASENAME_OUTER_EDGE_CYCLES);
}

static void dump_outer_edge_cycles_identified(const std::vector<outer_cycle>& outer_cycles,
                                              const std::vector<int>&         outer_edges,
                                              const v_CS&                     sequence_candidates,
                                              const std::vector<PointInt>&    points,

                                              const int* outer_cycle_pair, int iclockwise,
                                              const int* iedge_top)

{
    FILE* fp = fopen(DUMP_BASENAME_IDENTIFIED_OUTER_EDGE_CYCLE, "w");
    assert(fp);

    // the kernel limits the #! line to 127 characters, so I abbreviate
    fprintf(fp, "#!/usr/bin/feedgnuplot --datai --dom --aut --square --rangesizea 3 --w 'vec size screen 0.01,20 fixed fill' --set 'yr [:] rev'\n");
    fprintf(fp, "# fromx type fromy deltax deltay\n");

    for( int i_cycle=0; i_cycle<2; i_cycle++ )
    {
        for(int i_edge = 0; i_edge<4; i_edge++ )
        {
            const CandidateSequence* cs = &sequence_candidates[outer_edges[ outer_cycles[outer_cycle_pair[i_cycle]].e[i_edge] ]];
            const PointInt*          pt = &points[cs->c0->source_index()];

            char what[128];
            sprintf(what, "%s%s",
                    (i_cycle == iclockwise) ? "clockwise" : "counterclockwise",
                    iedge_top[i_cycle] == i_edge ? "-top" : "");

            fprintf(fp,
                    "%f %s %f %f %f\n",
                    (double)(pt->x) / (double)FIND_GRID_SCALE,
                    what,
                    (double)(pt->y) / (double)FIND_GRID_SCALE,
                    cs->delta_mean.x / (double)FIND_GRID_SCALE,
                    cs->delta_mean.y / (double)FIND_GRID_SCALE);
        }
    }

    fclose(fp);
    chmod(DUMP_BASENAME_IDENTIFIED_OUTER_EDGE_CYCLE,
          S_IRUSR | S_IRGRP | S_IROTH |
          S_IWUSR | S_IWGRP |
          S_IXUSR | S_IXGRP | S_IXOTH);
    fprintf(stderr, "Wrote outer edge cycle dump to %s\n",
            DUMP_BASENAME_IDENTIFIED_OUTER_EDGE_CYCLE);
}

static bool is_crossing( unsigned int line0_pt0, unsigned int line0_pt1,
                         unsigned int line1_pt0, unsigned int line1_pt1,
                         const std::vector<PointInt>&    points)
{
    // First I translate everything so that line0_pt0 is at the origin
    float l0pt1[2] =
        { (float)(points[line0_pt1].x - points[line0_pt0].x),
          (float)(points[line0_pt1].y - points[line0_pt0].y) };
    float l1pt0[2] =
        { (float)(points[line1_pt0].x - points[line0_pt0].x),
          (float)(points[line1_pt0].y - points[line0_pt0].y) };
    float l1pt1[2] =
        { (float)(points[line1_pt1].x - points[line0_pt0].x),
          (float)(points[line1_pt1].y - points[line0_pt0].y) };

    // Now I rotate the 3 points such that l0pt1 ends up aligned on the x axis.
    // I don't bother to divide by the hypotenuse, so l0pt1 is at (d^2, 0)
    float d2 = l0pt1[0]*l0pt1[0] + l0pt1[1]*l0pt1[1];

    float l1pt0_rotated[2] =
        {  l1pt0[0]*l0pt1[0] + l1pt0[1]*l0pt1[1],
          -l1pt0[0]*l0pt1[1] + l1pt0[1]*l0pt1[0] };
    float l1pt1_rotated[2] =
        {  l1pt1[0]*l0pt1[0] + l1pt1[1]*l0pt1[1],
          -l1pt1[0]*l0pt1[1] + l1pt1[1]*l0pt1[0] };

    // In this rotated space, the two points must have opposite-sign y (lie on
    // both sides of the line).
    if( l1pt0_rotated[1]*l1pt1_rotated[1] > 0 )
        return false;

    // To cross the line, both x coords can't be off the edge
    if( (l1pt0_rotated[0] < 0  && l1pt1_rotated[0] < 0) ||
        (l1pt0_rotated[0] > d2 && l1pt1_rotated[0] > d2) )
        return false;

    // OK, maybe they do cross. I actually compute the intersection
    // a + k(b-a) = [..., 0] -> k = a1/(a1-b1)
    float k = l1pt0_rotated[1] / (l1pt0_rotated[1] - l1pt1_rotated[1]);
    float x = l1pt0_rotated[0] + k * (l1pt1_rotated[0] - l1pt0_rotated[0]);
    return x >= 0.0f && x <= d2;
}

// recursively search for 4-cycle sequences of outer edges. On success, returns
// true, with the cycle reported in edges
static bool next_outer_edge(// this iteration
                            outer_cycle* edges,      // edge sequence being
                            // evaluated. Output
                            // returned here
                            int          edge_count, // how many edges we've got,
                            // including this one

                            // context
                            unsigned int                    point_initial,
                            const std::vector<int>&         outer_edges,
                            const v_CS&                     sequence_candidates,
                            const SequenceIndicesFromPoint& outer_edges_from_point,
                            const std::vector<PointInt>&    points,
                            bool                            debug)
{
    bool        found_cycle = false;
    outer_cycle outer_cycle_found = {};

    int i_edge = (int)edges->e[edge_count-1];
    unsigned int first_point_this_edge = sequence_candidates[outer_edges[i_edge]].c0   ->source_index();
    unsigned int last_point_this_edge  = sequence_candidates[outer_edges[i_edge]].clast->source_index();

    const std::vector<int>* next_edges;
    try
    {
        next_edges = &outer_edges_from_point.at(last_point_this_edge);
    }
    catch(...)
    {
        if(debug)
            fprintf(stderr, "No opposing outer edge\n");
        return false;
    }
    int Nedges_from_here = next_edges->size();
    for(int i=0; i<Nedges_from_here; i++)
    {
        // I make sure to not follow edges that are inverses of the immediately
        // previous edge, and I make sure that the 4th edge forms a loop AND
        // that a non-4th edge does NOT go back to the start. This is enough,
        // and I don't need any more already-visited logic:
        //
        // - 1st edge is given
        // - 2nd edge can go anywhere except directly backwards
        // - 3rd edge can go anywhere except the edges 1 (the start) and 2 (the
        //   previous point)
        // - 4th edge can go only to the start point
        unsigned int last_point_next_edge = sequence_candidates[outer_edges[ (*next_edges)[i] ]].clast->source_index();

        if( last_point_next_edge == first_point_this_edge )
            // This next edge is an inverse of this edge. It's not a part of my
            // 4-cycle
            continue;

        // I need to ignore X structures. In the below, ACBDA is valid,
        // but ABCDA is NOT valid.
        //
        //   A  C
        //   |\/|
        //   |/\|
        //   D  B
        //
        // Thus I make sure that
        // edge 3 does not cross edge 1 and that
        // edge 4 does not cross edge 2
        if(edge_count != 3)
        {
            // This is not the last edge. It may not go to the start
            if( last_point_next_edge == point_initial )
                continue;

            if(edge_count == 2)
            {
                if( is_crossing(sequence_candidates[outer_edges[ edges->e[0]      ]].c0   ->source_index(),
                                sequence_candidates[outer_edges[ edges->e[0]      ]].clast->source_index(),
                                sequence_candidates[outer_edges[ (*next_edges)[i] ]].c0   ->source_index(),
                                sequence_candidates[outer_edges[ (*next_edges)[i] ]].clast->source_index(),
                                points ))
                    continue;
            }

            edges->e[edge_count] = (short)(*next_edges)[i];
            if(!next_outer_edge( edges, edge_count+1,
                                 point_initial,
                                 outer_edges,
                                 sequence_candidates,
                                 outer_edges_from_point,
                                 points,
                                 debug))
                continue;

            // Got a successful result. It must be unique
            if(found_cycle)
            {
                if(debug)
                    fprintf(stderr, "Found non-unique 4-cycle\n");
                return false;
            }
            found_cycle = true;
            outer_cycle_found = *edges;
        }
        else
        {
            // Last edge. May only go back to the start
            if( last_point_next_edge != point_initial )
                continue;

            if( is_crossing(sequence_candidates[outer_edges[ edges->e[1]      ]].c0   ->source_index(),
                            sequence_candidates[outer_edges[ edges->e[1]      ]].clast->source_index(),
                            sequence_candidates[outer_edges[ (*next_edges)[i] ]].c0   ->source_index(),
                            sequence_candidates[outer_edges[ (*next_edges)[i] ]].clast->source_index(),
                            points ))
            {
                // I already found the last edge, but it's crossing itself. I
                // know there aren't any more solutions here, so I exit
                return false;
            }

            edges->e[3] = (short)(*next_edges)[i];
            return true;
        }
    }

    if(!found_cycle) return false;

    *edges = outer_cycle_found;
    return true;
}

static bool is_equalAndOpposite_cycle(const outer_cycle& cycle0,
                                      const outer_cycle& cycle1,

                                      // context
                                      const std::vector<int>& outer_edges,
                                      const v_CS&             sequence_candidates,
                                      bool                    debug)
{
    // Pick an arbitrary starting point: initial point of the first edge of the
    // first cycle
    int          iedge0 = 0;
    unsigned int ipt0   = sequence_candidates[outer_edges[cycle0.e[iedge0]]].c0->source_index();

    // find the edge in the potentially-opposite cycle that ends at this point
    int iedge1 = -1;
    for(int _iedge1 = 0; _iedge1 < 4; _iedge1++)
        if(ipt0 == sequence_candidates[outer_edges[ cycle1.e[_iedge1] ]].clast->source_index())
        {
            iedge1 = _iedge1;
            break;
        }
    if( iedge1 < 0)
    {
        if(debug)
            fprintf(stderr, "Given outer cycles are NOT equal and opposite: couldn't find a corresponding point in the two cycles\n");
        return false;
    }

    // Now I traverse the two cycles, and compare. I traverse the second cycle
    // backwards, and compare the back/front and front/back
    for(int i=0; i<4; i++)
    {
        unsigned int cycle0_points[2] =
            { (unsigned int)sequence_candidates[outer_edges[cycle0.e[iedge0]]].c0   ->source_index(),
              (unsigned int)sequence_candidates[outer_edges[cycle0.e[iedge0]]].clast->source_index()};
        unsigned int cycle1_points[2] =
            { (unsigned int)sequence_candidates[outer_edges[cycle1.e[iedge1]]].c0   ->source_index(),
              (unsigned int)sequence_candidates[outer_edges[cycle1.e[iedge1]]].clast->source_index() };
        if(cycle0_points[0] != cycle1_points[1] ||
           cycle0_points[1] != cycle1_points[0] )
        {
            if(debug)
                fprintf(stderr, "Given outer cycles are NOT equal and opposite\n");
            return false;
        }

        iedge0 = (iedge0+1) % 4;
        iedge1 = (iedge1+3) % 4;
    }
    return true;
}

// assumes the two given cycles are equal and opposite
//
// returns
//   0  if cycle0 is clockwise
//   1  if cycle1 is clockwise
//   <0 if the cycle is not convex
//
// Convexity is determined by connecting the corners with straight lines. I
// guess maybe a lens could be so distorted that you get nonconvex polygon, but
// I can't quite imagine that. The polygon would have to look like this:
/*

            /\
           /  \
          /    \
         /      \
        /   /\   \
       /  --  --  \
      /__/      \__\
*/
static int select_clockwise_cycle_and_find_top(// out
                                               int iedge_top[2],

                                               const outer_cycle& cycle0,
                                               const outer_cycle& cycle1,

                                               // context
                                               const std::vector<int>&      outer_edges,
                                               const v_CS&                  sequence_candidates,
                                               const std::vector<PointInt>& points,
                                               bool                         debug)
{
    // I pick a cycle, and look at the sign of the cross-product of each
    // consecutive direction. The sign should be constant for all consecutive
    // directions. Non-convexity would generate a different sign. And the sign
    // tells me if the thing is clockwise or not

    // 4 segments, each direction has an x and a y
    int v[4][2];
    for(int i=0; i<4; i++)
    {
        unsigned int ipt0 = sequence_candidates[outer_edges[cycle0.e[i]]].c0   ->source_index();
        unsigned int ipt1 = sequence_candidates[outer_edges[cycle0.e[i]]].clast->source_index();

        v[i][0] = (points[ipt1].x - points[ipt0].x) / FIND_GRID_SCALE_APPROX_POWER2 ;
        v[i][1] = (points[ipt1].y - points[ipt0].y) / FIND_GRID_SCALE_APPROX_POWER2 ;
    }

    bool sign[4];
    for(int i0=0; i0<4; i0++)
    {
        int i1 = (i0+1)%4;
        sign[i0] = v[i1][0]*v[i0][1] < v[i0][0]*v[i1][1];
    }

    int i_clockwise;
    if(      sign[0] &&  sign[1] &&  sign[2] &&  sign[3]) i_clockwise = 0;
    else if(!sign[0] && !sign[1] && !sign[2] && !sign[3]) i_clockwise = 1;
    else
    {
        if(debug)
            fprintf(stderr, "The outer edge cycles aren't convex!\n");
        return -1;
    }


    // To find the "top", I pick the edge with the lowest y coord of its center
    // point. Everything is assumed to be sorta square
    const outer_cycle* cycles[2] = {&cycle0, &cycle1};
    for( int icycle=0; icycle<2; icycle++)
    {
        int ymid2_min = INT_MAX;
        int iedge_min = -1;

        for(int i=0; i<4; i++)
        {
            unsigned int ipt0 = sequence_candidates[outer_edges[cycles[icycle]->e[i]]].c0   ->source_index();
            unsigned int ipt1 = sequence_candidates[outer_edges[cycles[icycle]->e[i]]].clast->source_index();

            int ymid2 = (int)(points[ipt0].y + points[ipt1].y);
            if(ymid2 < ymid2_min)
            {
                ymid2_min = ymid2;
                iedge_min = i;
            }
        }

        iedge_top[icycle] = iedge_min;
    }

    return i_clockwise;
}

static int find_sequence_from_to( // inputs
                                  unsigned int from, unsigned int to,

                                  // context
                                  const v_CS&  sequence_candidates,
                                  const SequenceIndicesFromPoint& sequences_from_point )
{
    try
    {
        const std::vector<int>& sequences = sequences_from_point.at(from);
        for(int i=0; i<(int)sequences.size(); i++)
        {
            const CandidateSequence* cs = &sequence_candidates[sequences[i]];
            if(cs->clast->source_index() == to)
                return sequences[i];
        }
        return -1;
    }
    catch(...)
    {
        return -1;
    }
}

__attribute__((visibility("default")))
bool mrgingham::find_grid_from_points( // out
                                      std::vector<PointDouble>& points_out,

                                      // in
                                      const std::vector<PointInt>& points,
                                      bool     debug,
                                      const debug_sequence_t& debug_sequence)
{
    VORONOI voronoi;
    construct_voronoi(points.begin(), points.end(), &voronoi);

    if(debug)
        dump_voronoi(&voronoi, points);

    v_CS sequence_candidates;
    get_sequence_candidates(&sequence_candidates, &voronoi, points,
                            debug_sequence);

    if(debug)
    {
        dump_candidates(DUMP_BASENAME_ALL_SEQUENCE_CANDIDATES,
                        &sequence_candidates, NULL, points);

        fprintf(stderr, "got %zd points\n", points.size());
        fprintf(stderr, "got %zd sequence candidates\n", sequence_candidates.size());
    }

    // I have all the sequence candidates. I find all the sequences that could
    // be edges of my grid: each one begins at a cell that's the start of at
    // least two sequences
    std::vector<int> outer_edges;
    // I likely only need 8, but I don't want to ever reallocate this thing
    outer_edges.reserve(20);
    std::map<unsigned int, int> sequences_initiated_count;
    int Ncs = sequence_candidates.size();
    for( int i=0; i<Ncs; i++ )
    {
        const CandidateSequence* cs = &sequence_candidates[i];
        sequences_initiated_count[cs->c0->source_index()]++;
    }
    for( int i=0; i<Ncs; i++ )
    {
        const CandidateSequence* cs = &sequence_candidates[i];
        if(sequences_initiated_count[cs->c0->source_index()] >= 2)
            outer_edges.push_back(i);
    }

    // I now have potential outer edges: all sequences that begin with a cell
    // that initiates at least two sequences (this one and at least one other)
    if( outer_edges.size() < 8 )
    {
        // need at least 8 outer edges: 4 in each direction
        if(debug)
        {
            fprintf(stderr, "Too few candidates for an outer edge of the grid. Needed at least 8, got %d\n",
                    (int)outer_edges.size());
        }
        return false;
    }
    if(debug)
        dump_candidates(DUMP_BASENAME_OUTER_EDGES,
                        &sequence_candidates, &outer_edges, points);

    // I won't have very many of these outer edges, so I don't worry too much
    // about efficient algorithms here.
    //
    // I look for 2 cyclical sequences of length 4 each. Equal and opposite to
    // each other
    int Nouter_edges = outer_edges.size();

    SequenceIndicesFromPoint outer_edges_from_point;
    for( int i=0; i<Nouter_edges; i++ )
    {
        const CandidateSequence* cs = &sequence_candidates[outer_edges[i]];
        outer_edges_from_point[cs->c0->source_index()].push_back(i);
    }

    std::vector<outer_cycle> outer_cycles;
    std::set<int> outer_edges_in_found_cycles;
    for( int i=0; i<Nouter_edges; i++ )
    {
        if( outer_edges_in_found_cycles.count(i) )
            // I already processed this edge
            continue;

        outer_cycle outer_cycle_found;
        outer_cycle_found.e[0] = i;
        if(!next_outer_edge(// this iteration
                            &outer_cycle_found,
                            1,

                            // context
                            sequence_candidates[outer_edges[i]].c0->source_index(),
                            outer_edges,
                            sequence_candidates,
                            outer_edges_from_point,
                            points,
                            debug))
            continue;

        outer_cycles.push_back( outer_cycle_found );
        for(int i=0; i<4; i++)
            outer_edges_in_found_cycles.insert(outer_cycle_found.e[i]);
    }

    if(debug && outer_cycles.size())
        dump_outer_edge_cycles(outer_cycles, outer_edges, sequence_candidates, points);

    if(outer_cycles.size() < 2)
    {
        if(debug)
            fprintf(stderr, "Found too few 4-cycles. Needed at least 2, got %d\n",
                    (int)outer_cycles.size());
        return false;
    }

    // I should have exactly one set of an equal/opposite cycles
    int outer_cycle_pair[2] = {-1,-1};
    for(int i0=0; i0<(int)outer_cycles.size(); i0++)
        for(int i1=i0+1; i1<(int)outer_cycles.size(); i1++)
        {
            if(is_equalAndOpposite_cycle(outer_cycles[i0], outer_cycles[i1],
                                         outer_edges, sequence_candidates,
                                         debug))
            {
                if(outer_cycle_pair[0] >= 0)
                {
                    if(debug)
                        fprintf(stderr, "Found more than one equal-and-opposite pair of outer-edge cycles. Giving up\n");
                    return false;
                }
                outer_cycle_pair[0] = i0;
                outer_cycle_pair[1] = i1;
            }
        }
    if(outer_cycle_pair[0] < 0)
    {
        if(debug)
            fprintf(stderr, "Didn't find any equal-and-opposite pairs of outer-edge cycles. Giving up\n");
        return false;
    }

    // I have my equal-and-opposite pair of cycles. I find the clockwise one. It
    // contains the top edge, which is the first in the sequence I'm going to
    // end up reporting
    int iedge_top[2];
    int iclockwise =
        select_clockwise_cycle_and_find_top(// out
                                            iedge_top,

                                            // cycles I'm looking at
                                            outer_cycles[outer_cycle_pair[0]],
                                            outer_cycles[outer_cycle_pair[1]],

                                            // context
                                            outer_edges, sequence_candidates,
                                            points,
                                            debug);
    if(iclockwise < 0)
        return false;

    if(debug)
        dump_outer_edge_cycles_identified(outer_cycles, outer_edges, sequence_candidates,
                                          points,
                                          outer_cycle_pair, iclockwise,
                                          iedge_top);

    // All done with the outer edges of the board. I now fill-in the internal
    // grid
    SequenceIndicesFromPoint sequences_from_point;
    for( int i=0; i<(int)sequence_candidates.size(); i++ )
    {
        const CandidateSequence* cs = &sequence_candidates[i];
        sequences_from_point[cs->c0->source_index()].push_back(i);
    }

    // sequences in sequence_candidates[]
    int horizontal_rows[Nwant];
    int vertical_left, vertical_right;

    horizontal_rows[0] = outer_edges[outer_cycles[outer_cycle_pair[  iclockwise]].e[  iedge_top[  iclockwise]          ]];
    vertical_left      = outer_edges[outer_cycles[outer_cycle_pair[1-iclockwise]].e[ (iedge_top[1-iclockwise] + 1) % 4 ]];
    vertical_right     = outer_edges[outer_cycles[outer_cycle_pair[  iclockwise]].e[ (iedge_top[  iclockwise] + 1) % 4 ]];

    unsigned int vertical_left_points [Nwant];
    unsigned int vertical_right_points[Nwant];
    get_candidate_points( vertical_left_points,  &sequence_candidates[vertical_left ], points );
    get_candidate_points( vertical_right_points, &sequence_candidates[vertical_right], points );

    for(int i=1; i<Nwant; i++)
    {
        // I fill in horizontal_rows[i]. I know each row must start at
        // vertical_left[i] and end at vertical_right[i]
        int sequence =
            find_sequence_from_to( vertical_left_points[i], vertical_right_points[i],
                                   sequence_candidates, sequences_from_point );

        if( sequence < 0 )
        {
            if(debug)
                fprintf(stderr, "Couldn't find sequence in row %d\n", i);
            return false;
        }

        horizontal_rows[i] = sequence;

        // Let's make sure the sequence from the other direction also works
        sequence =
            find_sequence_from_to( vertical_right_points[i], vertical_left_points[i],
                                   sequence_candidates, sequences_from_point );
        if(sequence < 0)
        {
            if(debug)
                fprintf(stderr, "Row %d: left-to-right sequence was found, but right-to-left sequence doesn't exist!\n", i);
            return false;
        }
    }

    // DO AGAIN AS A TRANSPOSED THING TO CONFIRM

    for(int i=0; i<Nwant; i++)
        output_row(points_out, sequence_candidates[horizontal_rows[i]], points);

    if(debug)
        fprintf(stderr, "Success. Found grid\n");
    return true;
}
