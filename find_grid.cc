#include <sys/stat.h>
#include <stdio.h>
#include <vector>
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



#define FOR_ALL_ADJACENT_CELLS(c) do {                                  \
    const PointInt*        pt = &points[c->source_index()];                \
    const VORONOI::edge_type* const e0 = c->incident_edge();            \
    bool first = true;                                                  \
    for(const VORONOI::edge_type* e = e0;                               \
        e0 != NULL && (e != e0 || first);                               \
        e = e->next(), first=false)                                     \
    {                                                                   \
        const VORONOI::cell_type* c_adjacent = e->twin()->cell();       \
                                                                        \
        const PointInt* pt_adjacent = &points[c_adjacent->source_index()]; \
                                                                        \
        PointInt delta( pt_adjacent->x - pt->x,                            \
                        pt_adjacent->y - pt->y );

#define FOR_ALL_ADJACENT_CELLS_END() }} while(0)






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
            fprintf(stderr, "..... accepting!\n");
        return c_adjacent;

    } FOR_ALL_ADJACENT_CELLS_END();

    return NULL;
}

static bool search_along_sequence( // out
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

    FOR_MATCHING_ADJACENT_CELLS(debug_sequence_pointscale)
    {
        if( c_adjacent == NULL )
            return false;
        delta_mean->x += (double)stats.delta_last.x;
        delta_mean->y += (double)stats.delta_last.y;
    }
    FOR_MATCHING_ADJACENT_CELLS_END();

    delta_mean->x /= (double)(N_remaining+1);
    delta_mean->y /= (double)(N_remaining+1);

    return true;
}

static void write_cell_center( std::vector<PointDouble>& points_out,
                               const VORONOI::cell_type* c,
                               const std::vector<PointInt>& points )
{
    const PointInt* pt = &points[c->source_index()];
    points_out.push_back( PointDouble( (double)pt->x / (double)FIND_GRID_SCALE,
                                       (double)pt->y / (double)FIND_GRID_SCALE) );
}

static void write_along_sequence( std::vector<PointDouble>& points_out,
                                  const PointInt* delta,
                                  const VORONOI::cell_type* c,
                                  int N_remaining,

                                  const std::vector<PointInt>& points)
{
    FOR_MATCHING_ADJACENT_CELLS(-1)
    {
        write_cell_center(points_out, c_adjacent, points);
    } FOR_MATCHING_ADJACENT_CELLS_END();
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
                                           int i_candidate,
                                           const PointInt* delta,
                                           const VORONOI::cell_type* c,
                                           int N_remaining,

                                           const std::vector<PointInt>& points)
{
    FOR_MATCHING_ADJACENT_CELLS(-1)
    {
        dump_interval(fp, i_candidate, i+1, c, c_adjacent, points);
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

struct CandidateSequence
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

typedef std::vector<CandidateSequence> v_CS;

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
                fprintf(stderr, "====== Looking at adjacent point (%d,%d)\n",
                        pt_adjacent->x / debug_sequence_pointscale,
                        pt_adjacent->y / debug_sequence_pointscale);

            PointDouble delta_mean;
            if( search_along_sequence( &delta_mean,
                                       &delta, c_adjacent, Nwant-2, points,
                                       (c == tracing_c) ? debug_sequence_pointscale : -1 ) )
            {
                double spacing_angle  = get_spacing_angle(delta_mean.y, delta_mean.x);
                double spacing_length = hypot(delta_mean.x, delta_mean.y);

                sequence_candidates->push_back( CandidateSequence({c, c_adjacent, delta_mean,
                                                                   spacing_angle, spacing_length}) );
            }
        } FOR_ALL_ADJACENT_CELLS_END();
    }
}

struct ClassificationBin
{
    PointDouble delta_mean_sum;
    int N;
};

#define THRESHOLD_BINFIT_LENGTH (120.0*(double)FIND_GRID_SCALE)
#define THRESHOLD_BINFIT_ANGLE  40.0

static bool fits_in_bin( const CandidateSequence* cs,
                         const ClassificationBin* bin,

                         double threshold_binfit_length,
                         double threshold_binfit_angle)
{
    if( bin->N == 0 ) return true;

    double dx         = bin->delta_mean_sum.x/(double)bin->N;
    double dy         = bin->delta_mean_sum.y/(double)bin->N;
    double bin_length = hypot(dx,dy);
    double bin_angle  = get_spacing_angle(dy, dx);

    double length_err = cs->spacing_length - bin_length;
    if(length_err*length_err > threshold_binfit_length*threshold_binfit_length) return false;

    double angle_err = cs->spacing_angle - bin_angle;

    // I want the angular error to be in [-90:90]
    angle_err += 10.0 * 180.0;
    angle_err = remainder(angle_err, 180.0);

    if(angle_err*angle_err > threshold_binfit_angle*threshold_binfit_angle) return false;

    return true;
}

static void push_to_bin( CandidateSequence* cs,
                         ClassificationBin* bin,
                         int bin_index )
{
    // I want to accumulate the vector with a consistent absolute direction
    // (angle modu0 180)
    if(bin->delta_mean_sum.x * cs->delta_mean.x +
       bin->delta_mean_sum.y * cs->delta_mean.y >= 0.0 )
    {
        bin->delta_mean_sum.x += cs->delta_mean.x;
        bin->delta_mean_sum.y += cs->delta_mean.y;
    }
    else
    {
        bin->delta_mean_sum.x -= cs->delta_mean.x;
        bin->delta_mean_sum.y -= cs->delta_mean.y;
    }

    cs->bin_index_neg = -bin_index - 1;
    bin->N++;
}

static int gather_unclassified(ClassificationBin* bin,
                               v_CS* sequence_candidates,
                               int bin_index)
{
    int Nremaining = 0;

    *bin = ClassificationBin({});

    for( auto it = sequence_candidates->begin(); it != sequence_candidates->end(); it++ )
    {
        CandidateSequence* cs = &(*it);

        if( cs->type != UNCLASSIFIED )
            continue;

        if( fits_in_bin(cs, bin, THRESHOLD_BINFIT_LENGTH, THRESHOLD_BINFIT_ANGLE) )
            push_to_bin(cs, bin, bin_index);
        else
            Nremaining++;
    }

    return Nremaining;
}

static void mark_outliers( v_CS* sequence_candidates,
                           int bin_index )
{
    for( auto it = sequence_candidates->begin(); it != sequence_candidates->end(); it++ )
    {
        CandidateSequence* cs = &(*it);

        if( // if we're calling everything remaining and outlier
           (bin_index < 0 && cs->type == UNCLASSIFIED) ||

            // or we're looking at THIS bin
           (cs->bin_index_neg == -bin_index - 1) )
        {
            cs->type = OUTLIER;
        }
    }
}

static void mark_orientation( v_CS* sequence_candidates,
                              const enum ClassificationType* types )
{
    for( auto it = sequence_candidates->begin(); it != sequence_candidates->end(); it++ )
    {
        CandidateSequence* cs = &(*it);

        if(      cs->bin_index_neg == -1 ) cs->type = types[0];
        else if( cs->bin_index_neg == -2 ) cs->type = types[1];
    }
}

static bool cluster_sequence_candidates( v_CS* sequence_candidates )
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
    ClassificationBin bins[3];

    int bin_index = 0;
    while(true)
    {
        ClassificationBin* bin = &bins[bin_index];
        int Nremaining = gather_unclassified( bin, sequence_candidates, bin_index );

        if( bin->N < Nwant*2 ) // should have enough for both directions
        {
            // this is a bin full of outliers
            mark_outliers( sequence_candidates, bin_index );
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
            mark_outliers(sequence_candidates, -1);
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

    mark_orientation( sequence_candidates, bin_orientation );
    return true;
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

static bool compare_reverse_along_sequence( const unsigned int* cs_points_other,

                                            const PointInt* delta,
                                            const VORONOI::cell_type* c,
                                            int N_remaining,

                                            const std::vector<PointInt>& points)
{
    FOR_MATCHING_ADJACENT_CELLS(-1)
    {
        if(*cs_points_other != c_adjacent->source_index())
            return false;
        cs_points_other--;
    } FOR_MATCHING_ADJACENT_CELLS_END();

    return true;
}
static bool is_reverse_sequence( const unsigned int* cs_points_other,
                                 const CandidateSequence* cs,
                                 const std::vector<PointInt>& points )
{
    if( cs->c0->source_index() != cs_points_other[Nwant-1] ) return false;
    if( cs->c1->source_index() != cs_points_other[Nwant-2] ) return false;

    const PointInt* pt0 = &points[cs->c0->source_index()];
    const PointInt* pt1 = &points[cs->c1->source_index()];

    PointInt delta({ pt1->x - pt0->x,
                  pt1->y - pt0->y});
    return compare_reverse_along_sequence(&cs_points_other[Nwant-3], &delta, cs->c1, Nwant-2, points);
}

static bool matches_direction(CandidateSequence* cs,
                              ClassificationType orientation )
{
    if( orientation == HORIZONTAL ) return cs->delta_mean.x > 0.0;
    return cs->delta_mean.y > 0.0;
}

static void filter_bidirectional( v_CS* sequence_candidates,
                                  const std::vector<PointInt>& points,
                                  ClassificationType orientation )
{
    // I loop through the candidates list, and try to find a matching other
    // candidate that is THIS candidate in the opposite order.
    //
    // If no such match is found, I throw away the candidate.
    // If such match IS found, I throw away one of the two
    int N = sequence_candidates->size();
    for( int i=0; i<N; i++ )
    {
        CandidateSequence* cs0 = &(*sequence_candidates)[i];
        if(cs0->type != orientation) continue;

        unsigned int cs0_points[Nwant];
        get_candidate_points(cs0_points, cs0, points);

        bool found = false;
        for( int j=i+1; j<N; j++ )
        {
            CandidateSequence* cs1 = &(*sequence_candidates)[j];
            if(cs1->type != orientation) continue;

            if( !is_reverse_sequence( cs0_points, cs1, points ) )
                continue;

            // bam. found reverse sequence. Throw away one of the matches. I
            // keep the one that matches the canonical direction the best ([1,0]
            // for "horizontal" and [0,1] for "vertical")
            if( !matches_direction(cs0, orientation) )
                *cs0 = *cs1;
            cs1->type = OUTLIER;
            found = true;
            break;
        }
        if( !found )
            // this candidate doesn't have a match. Throw out self
            cs0->type = OUTLIER;
    }
}


// Looks through our classification and determines whether things look valid or
// not. Makes no changes to anything
static bool validate_clasification(const v_CS* sequence_candidates)
{
    // I should have exactly Nwant horizontal lines and Nwant vertical lines
    int Nhorizontal = 0;
    int Nvertical   = 0;

    for( auto it = sequence_candidates->begin(); it != sequence_candidates->end(); it++ )
    {
        if(      it->type == HORIZONTAL ) Nhorizontal++;
        else if( it->type == VERTICAL   ) Nvertical++;
    }
    if( Nhorizontal != Nwant ) return false;
    if( Nvertical   != Nwant ) return false;


    // OK then. The horizontal lines should each

    // I'm tired. Let's call this good enough for now
#warning complete
    return true;

}

// dumps a terse self-plotting vnlog visualization of sequence candidates, and a
// more detailed vnlog containing more data
#define DUMP_FILENAME_SEQUENCE_CANDIDATES_SPARSE_BEFORE "/tmp/mrgingham-3-candidates.vnl"
#define DUMP_FILENAME_SEQUENCE_CANDIDATES_DENSE_BEFORE  "/tmp/mrgingham-3-candidates-detailed.vnl"
#define DUMP_FILENAME_SEQUENCE_CANDIDATES_SPARSE_AFTER  "/tmp/mrgingham-4-candidates.vnl"
#define DUMP_FILENAME_SEQUENCE_CANDIDATES_DENSE_AFTER   "/tmp/mrgingham-4-candidates-detailed.vnl"
static void dump_candidates(const v_CS* sequence_candidates,
                            const std::vector<PointInt>& points,
                            bool post_filter)
{
    const char* dump_filename_sequence_candidates_sparse = post_filter ?
        DUMP_FILENAME_SEQUENCE_CANDIDATES_SPARSE_AFTER :
        DUMP_FILENAME_SEQUENCE_CANDIDATES_SPARSE_BEFORE;
    const char* dump_filename_sequence_candidates_dense = post_filter ?
        DUMP_FILENAME_SEQUENCE_CANDIDATES_DENSE_AFTER :
        DUMP_FILENAME_SEQUENCE_CANDIDATES_DENSE_BEFORE;

    FILE* fp = fopen(dump_filename_sequence_candidates_sparse, "w");
    assert(fp);

    // the kernel limits the #! line to 127 characters, so I abbreviate
    fprintf(fp, "#!/usr/bin/feedgnuplot --datai --dom --aut --square --rangesizea 3 --w 'vec size screen 0.01,20 fixed fill' --set 'yr [:] rev'\n");
    fprintf(fp, "# fromx type fromy deltax deltay\n");

    for( auto it = sequence_candidates->begin(); it != sequence_candidates->end(); it++ )
    {
        const CandidateSequence* cs = &(*it);
        const PointInt*             pt = &points[cs->c0->source_index()];

        fprintf(fp,
                "%f %s %f %f %f\n",
                (double)(pt->x) / (double)FIND_GRID_SCALE,
                type_string(cs->type),
                (double)(pt->y) / (double)FIND_GRID_SCALE,
                cs->delta_mean.x / (double)FIND_GRID_SCALE,
                cs->delta_mean.y / (double)FIND_GRID_SCALE);
    }
    fclose(fp);
    chmod(dump_filename_sequence_candidates_sparse,
          S_IRUSR | S_IRGRP | S_IROTH |
          S_IWUSR | S_IWGRP |
          S_IXUSR | S_IXGRP | S_IXOTH);
    fprintf(stderr, "Wrote self-plotting sequence-candidate dump to %s\n",
            dump_filename_sequence_candidates_sparse);


    // detailed
    fp = fopen(dump_filename_sequence_candidates_dense, "w");
    assert(fp);

    fprintf(fp, "# candidateid pointid fromx fromy tox toy deltax deltay len angle\n");
    int N = sequence_candidates->size();
    for( int i=0; i<N; i++ )
    {
        const CandidateSequence* cs = &(*sequence_candidates)[i];

        dump_interval(fp, i, 0, cs->c0, cs->c1, points);

        const PointInt* pt0 = &points[cs->c0->source_index()];
        const PointInt* pt1 = &points[cs->c1->source_index()];

        PointInt delta({ pt1->x - pt0->x,
                      pt1->y - pt0->y});
        dump_intervals_along_sequence( fp, i, &delta, cs->c1, Nwant-2, points);
    }
    fclose(fp);
    fprintf(stderr, "Wrote detailed sequence-candidate dump to %s\n",
            dump_filename_sequence_candidates_dense);
}

static void write_output( std::vector<PointDouble>& points_out,
                          const v_CS* sequence_candidates,
                          const std::vector<PointInt>& points )
{
    for( auto it = sequence_candidates->begin(); it != sequence_candidates->end(); it++ )
    {
        if( it->type == HORIZONTAL )
        {
            write_cell_center(points_out, it->c0, points);
            write_cell_center(points_out, it->c1, points);

            const PointInt* pt0 = &points[it->c0->source_index()];
            const PointInt* pt1 = &points[it->c1->source_index()];

            PointInt delta({ pt1->x - pt0->x,
                          pt1->y - pt0->y});
            write_along_sequence( points_out, &delta, it->c1, Nwant-2, points);
        }
    }
}

static void sort_candidates(v_CS* sequence_candidates,
                            const std::vector<PointInt>& points )
{
    // I sort my vertical sequences in order of increasing x
    //
    // I sort my horizontal sequences in order of increasing y


    struct S{
        bool operator() ( const CandidateSequence& a,
                          const CandidateSequence& b) const
        {
            if( a.type != b.type )
            {
                // HORIZONTAL is 1st, VERTICAL is 2nd, and I don't care about the others
                if( a.type == HORIZONTAL ) return true;
                if( b.type == HORIZONTAL ) return false;
                if( a.type == VERTICAL   ) return true;
                if( b.type == VERTICAL   ) return false;
                return a.type < b.type;
            }

            if( a.type == HORIZONTAL )
                return _points[a.c0->source_index()].y < _points[b.c0->source_index()].y;
            return _points[a.c0->source_index()].x < _points[b.c0->source_index()].x;
        }

        const std::vector<PointInt>& _points;
        S(const std::vector<PointInt>& __points) : _points(__points) {}
    } sequence_comparator(points);

    std::sort( sequence_candidates->begin(), sequence_candidates->end(),
               sequence_comparator );
}

static CandidateSequence* get_first(v_CS* sequence_candidates,
                                    ClassificationType orientation)
{
    int N = sequence_candidates->size();
    for( int i=0; i<N; i++ )
    {
        CandidateSequence* cs = &(*sequence_candidates)[i];
        if( cs->type == orientation ) return cs;
    }
    return NULL;
}

static bool filter_bounds(v_CS* sequence_candidates,
                          ClassificationType orientation,
                          const std::vector<PointInt>& points)
{
    // I look at the first horizontal sequence and make sure that it consists of
    // the first points of all the vertical sequences, in order. And vice versa
    //
    // This function will mark extra sequences as outliers, and it will return
    // false if anything is missing
    ClassificationType orientation_other;
    if( orientation == HORIZONTAL ) orientation_other = VERTICAL;
    else                            orientation_other = HORIZONTAL;

    CandidateSequence* cs_ref    = get_first(sequence_candidates, orientation);
    CandidateSequence* cs_others = get_first(sequence_candidates, orientation_other);
    if( cs_ref    == NULL ) return false;
    if( cs_others == NULL ) return false;

    unsigned int cs_ref_points[Nwant];
    get_candidate_points( cs_ref_points, cs_ref, points );
    int i;
    for(i=0; i<Nwant; i++, cs_others++)
    {
        if( cs_others->type != orientation_other )
            // no more valid other sequences to follow
            break;

        if( cs_ref_points[i] != cs_others->c0->source_index() )
        {
            // mismatch! One of these sequences is an outlier
#warning handle this
            return false;
        }
    }
    return i == Nwant;
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
        dump_candidates(&sequence_candidates, points, false);

        fprintf(stderr, "got %zd points\n", points.size());
        fprintf(stderr, "got %zd sequence candidates\n", sequence_candidates.size());
    }

    if( !cluster_sequence_candidates(&sequence_candidates))
    {
        if(debug)
            fprintf(stderr, "cluster_sequence_candidates() failed. No grid detected\n");
        return false;
    }

    filter_bidirectional(&sequence_candidates, points, HORIZONTAL);
    filter_bidirectional(&sequence_candidates, points, VERTICAL);

    if(debug)
        dump_candidates(&sequence_candidates, points, true);

    // This is relatively slow (I'm moving lots of stuff around by value), but
    // I'm likely to not feel it anyway
    sort_candidates(&sequence_candidates, points);

    if( !filter_bounds(&sequence_candidates, HORIZONTAL, points) )
    {
        if(debug)
            fprintf(stderr, "Horizontal sequence candidates out of bounds. No grid detected\n");
        return false;
    }
    if( !filter_bounds(&sequence_candidates, VERTICAL,   points) )
    {
        if(debug)
            fprintf(stderr, "Vertical sequence candidates out of bounds. No grid detected\n");
        return false;
    }
    if(!validate_clasification(&sequence_candidates))
    {
        if(debug)
            fprintf(stderr, "validate_clasification() failed. No grid detected\n");
        return false;
    }

    write_output(points_out, &sequence_candidates, points);
    if(debug)
        fprintf(stderr, "Success. Found grid\n");
    return true;
}
