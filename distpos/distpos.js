var ACTOR_RADIUS = 10;
var MOVESPEED_MULT = 20;

var clone = function(x) { return JSON.parse(JSON.stringify(x)); }
var clip = function(x, lo, hi) { return Math.max(lo, Math.min(hi, x)); };

// Matrix functions.
var zeros = function(n, m) {
    var mat = [];
    for (var i = 0; i < n; ++i) mat.push(new Array(m).fill(0));
    return mat;
}
var eye = function(n) {
    var mat = zeros(n, n);
    for (var i = 0; i < n; ++i) mat[i][i] = 1.0;
    return mat;
}

var madd = function(a, b) {
    console.assert(a.length == b.length);
    console.assert(a[0].length == b[0].length);
    var r = zeros(a.length, a[0].length);
    for (var i = 0; i < a.length; i++)
        for (var j = 0; j < a[0].length; j++)
            r[i][j] = a[i][j] + b[i][j];
    return r;
}
var msub = function(a, b) {
    console.assert(a.length == b.length);
    console.assert(a[0].length == b[0].length);
    var r = zeros(a.length, a[0].length);
    for (var i = 0; i < a.length; i++)
        for (var j = 0; j < a[0].length; j++)
            r[i][j] = a[i][j] - b[i][j];
    return r;
}
var mmult = function(a, b) {
    console.assert(a[0].length == b.length);
    var r = zeros(a.length, b[0].length);
    for (var i = 0; i < a.length; i++)
        for (var j = 0; j < b[0].length; j++)
            for (var k = 0; k < a[0].length; k++)
                r[i][j] += a[i][k] * b[k][j];
    return r;
}
var madds = function(a, b) {
    var r = zeros(a.length, a[0].length);
    for (var i = 0; i < a.length; i++)
        for (var j = 0; j < a[0].length; j++)
            r[i][j] = a[i][j] + b;
    return r;
}
var mmults = function(a, b) {
    var r = zeros(a.length, a[0].length);
    for (var i = 0; i < a.length; i++)
        for (var j = 0; j < a[0].length; j++)
            r[i][j] = a[i][j] * b;
    return r;
}

var mT = function(m) {
    var r = zeros(m[0].length, m.length);
    for (var i = 0; i < m.length; i++)
        for (var j = 0; j < m[0].length; j++)
            r[j][i] = m[i][j];
    return r;
}

// Moore-Penrose inverse.
var pinv = function(m, zero_eps=1e-8) {
    var svd = SVDJS.SVD(m);
    for (var i = 0; i < m.length; ++i)
        if (svd.q[i] > zero_eps)
            svd.q[i] = 1 / svd.q[i];
    for (var i = 0; i < m.length; ++i)
        for (var j = 0; j < m[0].length; ++j)
            svd.v[i][j] *= svd.q[j];
    return mmult(svd.v, mT(svd.u));
}

var blend = function(x) { return x * x * (3 - 2*x); }

var animate = function(t, ease_in, ease_out) {
    if (t <= ease_in) {
        return blend(t / ease_in);
    }
    t -= ease_in;
    if (t <= ease_out) {
        return blend(1 - t / ease_out);
    }
    return 0;
}

var uniform_random = function(lo, hi) { return lo + Math.random() * (hi - lo); };
var range_random = function(start, stop) {
    return clip(Math.floor(uniform_random(start, stop)),
                start, stop - 1);
}
var normal_random = function() {
    var u = 1 - Math.random(); // Converting [0,1) to (0, 1].
    var v = Math.random();
    return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
}

var draw_circle = function(ctx, x, y, r, fill) {
    ctx.beginPath();
    ctx.arc(x, y, r, 0, 2 * Math.PI);
    if (fill) {
        ctx.fill();
    } else {
        ctx.stroke();
    }
};


var read_config = function(form, canvas) {
    var form = document.getElementById("distpos-controls");
    var data = new FormData(form);
    var mindeg = Math.max(1, parseInt(data.get('mindeg')));
    var maxdeg = Math.max(mindeg, parseInt(data.get('maxdeg')));
    return {
        w: canvas.width,
        h: canvas.height,
        n: Math.max(0, parseInt(data.get('numnodes'))),
        connectivity: data.get('connectivity'),
        distcalc: data.get('distcalc'),
        mindeg: mindeg,
        maxdeg: maxdeg,
        movespeed: parseFloat(data.get('movespeed')),
        obsspeed: Math.max(0.1, parseFloat(data.get('obsspeed'))),
        iterspeed: Math.max(0.1, parseFloat(data.get('iterspeed'))),
        distnoise: Math.max(0.0, parseFloat(data.get('distnoise'))),
        weightdecay: Math.min(Math.max(0.0, parseFloat(data.get('weightdecay'))), 1.0),
        learnrate: Math.min(Math.max(0.0, parseFloat(data.get('learnrate'))), 1.0),
        running: data.get('running') !== null,
        drawarrows: data.get('drawarrows') !== null,
    };
};

var generate_actors = function(cfg, cur_actors) {
    var actors = clone(cur_actors);
    while (actors.length > cfg.n) actors.pop();
    while (actors.length < cfg.n) {
        actors.push({
            x: uniform_random(ACTOR_RADIUS, cfg.w - ACTOR_RADIUS),
            y: uniform_random(ACTOR_RADIUS, cfg.h - ACTOR_RADIUS),
            vx: normal_random() * MOVESPEED_MULT,
            vy: normal_random() * MOVESPEED_MULT,
        });
    }

    return actors;
};

var update_actors = function(cfg, actors, dt) {
    for (var i = 0; i < actors.length; ++i) {
        actors[i].x += dt * cfg.movespeed * actors[i].vx;
        actors[i].y += dt * cfg.movespeed * actors[i].vy;
        if (actors[i].x < ACTOR_RADIUS ||
            actors[i].x > cfg.w - ACTOR_RADIUS) actors[i].vx = -actors[i].vx;
        if (actors[i].y < ACTOR_RADIUS ||
            actors[i].y > cfg.h - ACTOR_RADIUS) actors[i].vy = -actors[i].vy;
        actors[i].x = clip(actors[i].x, ACTOR_RADIUS, cfg.w - ACTOR_RADIUS);
        actors[i].y = clip(actors[i].y, ACTOR_RADIUS, cfg.h - ACTOR_RADIUS);
    }
}

var generate_edges = function(cfg, actors, cur_estimate) {
    var degrees = [];
    while (degrees.length < cfg.n) {
        degrees.push(range_random(cfg.mindeg, cfg.maxdeg + 1));
    }

    var edgelist = [];
    if (cfg.connectivity == 'random') {
        for (var i = 0; i < cfg.n; ++i) {
            var edges = [];
            var edge_exists = {};
            edge_exists[i] = true; // Prevent self-edges.

            for (var d = 0; d < degrees[i]; ++d) {
                var j;
                do {
                    j = range_random(0, cfg.n);
                } while (edge_exists[j] === true);
                edges.push(j);
                edge_exists[j] = true;
            }
            edgelist.push(edges);
        }
    } else {
        // This can be more efficient than computing all distance pairs.
        var distances = [];
        for (var i = 0; i < cfg.n; ++i) {
            edgelist.push([]);

            for (var j = i + 1; j < cfg.n; ++j) {
                var dx, dy;
                if (cfg.distcalc == 'real') {
                    dx = actors[i].x - actors[j].x;
                    dy = actors[i].y - actors[j].y;
                } else {
                    dx = cur_estimate[i][0] - cur_estimate[j][0];
                    dy = cur_estimate[i][1] - cur_estimate[j][1];
                }
                var dist = dx*dx + dy*dy;
                distances.push([i, j, dist]);
            }
        }

        if (cfg.connectivity == 'closest') {
            distances.sort(function(a, b) { return a[2] - b[2]; });
        } else {
            distances.sort(function(a, b) { return b[2] - a[2]; });
        }

        for (var k = 0; k < distances.length; ++k) {
            var i = distances[k][0];
            var j = distances[k][1];
            if (edgelist[i].length < degrees[i]) edgelist[i].push(j);
            if (edgelist[j].length < degrees[j]) edgelist[j].push(i);
        }
    }

    return edgelist;
};

var generate_estimate = function(cfg) {
    var estimate = [];
    for (var i = 0; i < cfg.n; ++i) {
        estimate.push([uniform_random(0, cfg.w), uniform_random(0, cfg.h)]);
    }
    return estimate;
}

var compute_dists = function(cfg, pos, edgelist) {
    var distlist = [];
    for (var i = 0; i < edgelist.length; ++i) {
        var dists = [];
        for (var k = 0; k < edgelist[i].length; ++k) {
            var j = edgelist[i][k];

            var dx = pos[j].x - pos[i].x;
            var dy = pos[j].y - pos[i].y;
            var dist = Math.sqrt(dx*dx + dy*dy);
            dist *= 1 + normal_random() * cfg.distnoise;
            dists.push(dist);
        }

        distlist.push(dists);
    }

    return distlist;
}

var compute_dhat_w = function(edgelist, distlist) {
    var zero_eps = 1e-9;
    var n = edgelist.length;

    // Transform into distance matrix.
    var dhat = zeros(n, n);
    var w = zeros(n, n);
    for (var i = 0; i < n; ++i) {
        for (var k = 0; k < edgelist[i].length; ++k) {
            var j = edgelist[i][k];
            var d = distlist[i][k];

            if (i == j) {
                console.log("WARNING, i == j in edgelist");
                continue;
            }

            // Normalize (possibly conflicting) observations using average.
            if (dhat[j][i] > zero_eps) {
                d = (d + dhat[j][i]) / 2;
            }

            dhat[i][j] = d;
            dhat[j][i] = d;
            w[i][j] = 1;
            w[j][i] = 1;
        }
    }
    return [dhat, w];
}

var mds_v = function(w) {
    // Compute V.
    var n = w.length;
    var V = zeros(n, n);
    for (var i = 0; i < n; ++i) {
        for (var j = i + 1; j < n; ++j) {
            V[i][i] += w[i][j];
            if (i != j) {
                V[j][j] += w[i][j];
                V[i][j] -= w[i][j];
                V[j][i] -= w[i][j];
            }
        }
    }

    // Compute pseudo-inverse.
    var Vinv = pinv(V);

    return [V, Vinv];
}

var mds_step = function(cur_estimate, dhat, w, V, Vinv) {
    // SMACOF algorithm.
    // Implemented from Multidimensional Scaling by Majorization: A Review.
    // And sample implementation in numpy from https://github.com/scikit-learn/scikit-learn/issues/6828#issuecomment-510399353.

    var zero_eps = 1e-9;
    var n = cur_estimate.length;
    var Y = cur_estimate;

    // Compute B(Y).
    var BY = zeros(n, n);
    for (var i = 0; i < n; ++i) {
        for (var j = i + 1; j < n; ++j) {
            var dx = cur_estimate[j][0] - cur_estimate[i][0];
            var dy = cur_estimate[j][1] - cur_estimate[i][1];
            var d = Math.sqrt(dx*dx + dy*dy);
            var b = d <= zero_eps ? 0 : dhat[i][j] / d;

            BY[i][i] += w[i][j]*b;
            if (i != j) {
                BY[j][j] += w[i][j]*b;
                BY[i][j] -= w[i][j]*b;
                BY[j][i] -= w[i][j]*b;
            }
        }
    }

    return mmult(Vinv, mmult(BY, Y));
}


var correct_rotate_translate = function(estimate, truth) {
    // Least-Squares Rigid Motion Using SVD.
    // By Olga Sorkine-Hornung and Michael Rabinovich.
    var n = estimate.length;
    var truthmat = [];
    var estmat = [];
    var tc = [0, 0];
    var ec = [0, 0];

    // Compute centroids.
    for (var i = 0; i < n; ++i) {
        estmat.push([estimate[i][0], estimate[i][1]]);
        truthmat.push([truth[i].x, truth[i].y]);
        tc[0] += truthmat[i][0];
        tc[1] += truthmat[i][1];
        ec[0] += estmat[i][0];
        ec[1] += estmat[i][1];
    }

    tc[0] /= n; tc[1] /= n;
    ec[0] /= n; ec[1] /= n;
    
    // Subtract centroids.
    for (var i = 0; i < n; ++i) {
        truthmat[i][0] -= tc[0];
        truthmat[i][1] -= tc[1];
        estmat[i][0] -= ec[0];
        estmat[i][1] -= ec[1];
    }

    // Magic sauce.
    var H = mmult(mT(estmat), truthmat);
    var svd = SVDJS.SVD(H);
    var uT = mT(svd.u);
    var R = mmult(svd.v, uT);

    // var det = R[0][0]*R[1][1] - R[0][1]*R[1][0];
    // if (det < 0) {
    //     uT = mmult([[1, 0], [0, -1]], uT);
    //     R = mmult(svd.v, uT);
    // }
    // var det2 = R[0][0]*R[1][1] - R[0][1]*R[1][0];
    // console.log(det, det2);

    var rec = mmult(R, mT([ec]));
    var rot_est = mmult(R, mT(estimate));
    var trans = [tc[0] - rec[0], tc[1] - rec[1]];
    for (var i = 0; i < n; ++i) {
        rot_est[0][i] += trans[0];
        rot_est[1][i] += trans[1];
    }

    return mT(rot_est);
};






var distpos_demo = function() {
    var form = document.getElementById("distpos-controls");
    var canvas = document.getElementById("distpos-canvas");
    var ctx = canvas.getContext("2d");

    var cfg = read_config(form, canvas);
    var actors = generate_actors(cfg, []);
    var cur_estimate = generate_estimate(cfg);
    var dhat_t = zeros(cfg.n, cfg.n);
    var w_t = eye(cfg.n);
    var arrows = [];
    var last_t = null;
    var t_since_last_obs = 1/cfg.obsspeed + 1e-8;
    var t_since_last_step = 0;

    var randomize_estimate = function() {
        cur_estimate = generate_estimate(cfg);
        dhat_t = zeros(cfg.n, cfg.n);
        w_t = eye(cfg.n);
        t_since_last_obs = 1/cfg.obsspeed + 1e-8;
    };

    var randomize_agents = function() {
        actors = generate_actors(cfg, []);
    };

    form.addEventListener('change', function() {
        var ncfg = read_config(form, canvas);
        var cardinality_changed = ncfg.n !== cfg.n;
        var new_graph_needed = ncfg.n !== cfg.n
                            || ncfg.connectivity !== cfg.connectivity
                            || ncfg.mindeg !== cfg.mindeg
                            || ncfg.maxdeg !== cfg.maxdeg;

        cfg = ncfg;

        if (cardinality_changed) {
            actors = generate_actors(cfg, actors);
            randomize_estimate();
            arrows = [];
        }
    });

    document.getElementById("distpos-randagents")
        .addEventListener("click", randomize_agents, false);
    document.getElementById("distpos-randestimate")
        .addEventListener("click", randomize_estimate, false);

    var V, Vinv;

    var anim = function(t) {
        t /= 1000;
        var dt = last_t !== null ? t - last_t : 0;
        last_t = t;

        if (!cfg.running) dt = 0;

        update_actors(cfg, actors, dt);

        t_since_last_obs += dt;
        if (t_since_last_obs > 1/cfg.obsspeed) {
            t_since_last_obs -= 1/cfg.obsspeed;
            edgelist = generate_edges(cfg, actors, cur_estimate);
            var distlist = compute_dists(cfg, actors, edgelist);

            for (var i = 0; i < cfg.n; ++i) {
                for (var k = 0; k < edgelist[i].length; ++k) {
                    var j = edgelist[i][k];
                    arrows.push({t: t, start: i, stop: j});
                }
            }

            var [dhat, w] = compute_dhat_w(edgelist, distlist);
            for (var i = 0; i < cfg.n; ++i) {
                for (var j = 0; j < cfg.n; ++j) {
                    w_t[i][j] *= cfg.weightdecay;

                    if (w[i][j]) {
                        w_t[i][j] = 1;
                        dhat_t[i][j] = cfg.learnrate * dhat[i][j] + 
                                  (1 - cfg.learnrate) * dhat_t[i][j];
                    }
                }
            }

            [V, Vinv] = mds_v(w_t);
        }
        
        t_since_last_step += dt;
        while (t_since_last_step > 1/cfg.iterspeed) {
            t_since_last_step -= 1/cfg.iterspeed;
            cur_estimate = mds_step(cur_estimate, dhat_t, w_t, V, Vinv);
        }


        // Rendering
        ctx.clearRect(0, 0, cfg.w, cfg.h);

        ctx.font = "8px arial";
        ctx.textAlign = "center";
        ctx.textBaseline = "top";


        if (cfg.drawarrows) {
            // var alpha = 0.2 + 0.1*Math.pow(1 - t_since_last_obs / (1/cfg.obsspeed), 2);
            // var alpha = 0.5 * animate(t_since_last_obs, 0.05, 1 / (1/cfg.obsspeed) - 0.05);
            // var alpha = 0.2 * animate(t_since_last_obs, 0.001, 1/cfg.obsspeed - 0.001);
            for (var i = 0; i < arrows.length; ++i) {
                var arrow = arrows[i];
                var start = actors[arrow.start];
                var stop = actors[arrow.stop];
                var dx = stop.x - start.x;
                var dy = stop.y - start.y;
                var dist = Math.sqrt(dx*dx + dy*dy);

                var alpha = 1 / Math.cbrt(arrows.length) * animate(t - arrow.t, 0.1, 0.9);
                ctx.strokeStyle = 'rgb(0, 0, 0, ' + alpha + ')';
                ctx.fillStyle = ctx.strokeStyle;

                ctx.save();
                ctx.translate(start.x, start.y);
                ctx.rotate(Math.atan2(dy, dx));
                ctx.beginPath();
                ctx.moveTo(ACTOR_RADIUS, 0);

                var offset = ACTOR_RADIUS;
                ctx.lineTo(dist - offset - 10, 0);
                ctx.stroke();
                ctx.beginPath();
                ctx.lineTo(dist - offset - 10, -5);
                ctx.lineTo(dist - offset, 0);
                ctx.lineTo(dist - offset - 10, +5);

                ctx.fill();
                ctx.restore();
            }

            arrows = arrows.filter(function(arrow) {
                return t - arrow.t < 1;
            });
        }

        ctx.strokeStyle = 'black';
        ctx.fillStyle = 'black';
        ctx.textBaseline = 'top';
        for (var i = 0; i < actors.length; ++i) {
            ctx.fillStyle = 'white';
            draw_circle(ctx, actors[i].x, actors[i].y, ACTOR_RADIUS, true);
            ctx.fillStyle = 'black';
            draw_circle(ctx, actors[i].x, actors[i].y, ACTOR_RADIUS, false);
            ctx.fillText(1 + i, actors[i].x, actors[i].y + ACTOR_RADIUS + 3);
        }


        var display_estimate = correct_rotate_translate(cur_estimate, actors);
        ctx.strokeStyle = 'black';
        ctx.fillStyle = 'black';
        ctx.textBaseline = 'bottom';
        for (var i = 0; i < display_estimate.length; ++i) {
            draw_circle(ctx, display_estimate[i][0], display_estimate[i][1], 3, true);
            ctx.fillText(1 + i, display_estimate[i][0], display_estimate[i][1] - ACTOR_RADIUS - 3);
        }

        window.requestAnimationFrame(anim);
    }

    window.requestAnimationFrame(anim);
};

document.addEventListener('DOMContentLoaded', function() {
    distpos_demo();
}, false);
