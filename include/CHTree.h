//
// Created by etoga on 5/12/23.
//

#ifndef DYNAMICCONVEXHULL_CHTREE_H
#define DYNAMICCONVEXHULL_CHTREE_H

#include <vector>
#include <CGAL/enum.h>
#include "AvlTree.h"
#include "util.h"
#define isLeaf this->isLeaf


template<class Traits>
class CHTree : AVLTree<Bridges<Traits>>{
using Bridge = typename Traits::Segment_2;
using Point = typename Traits::Point_2;
using Node = typename AVLTree<Bridges<Traits>>::Node;
using Midpoint = typename Traits::Construct_midpoint_2;
using Compare_slope = typename Traits::Compare_slope_2;
using Compare_at_x = typename Traits::Compare_y_at_x_2;

const Midpoint midpoint = Midpoint();
const Compare_slope compare_slope = Compare_slope();
const Compare_at_x compare_at_x = Compare_at_x();

protected:
    template<bool lower>
    bool slope_comp(const Bridge& l, const Bridge& r){
        CGAL::Comparison_result res = compare_slope(l,r);
        if(res == CGAL::EQUAL) return true;
        return (res == CGAL::SMALLER) != lower;
    }

    template<bool lower>
    bool m_comp(const Bridge& l, const Bridge& r, const Point& m){
        if(l.is_vertical()) return !lower;
        if(r.is_vertical()) return lower;
        CGAL::Comparison_result res = compare_at_x(m,l.supporting_line(),r.supporting_line());
        if(res == CGAL::EQUAL) return true;
        return (res == CGAL::SMALLER) == lower;
    }

    template<bool lower>
    bool cover_comp(const Point& p, const Bridge& b){
        CGAL::Comparison_result res = compare_at_x(p,b.supporting_line());
        if(res == CGAL::EQUAL) return true;
        return (res == CGAL::SMALLER) == lower;
    }

    template<bool lower>
    Bridge findBridge(Node* v){
        Node* x = v->left;
        Node* y = v->right;
        Bridge e_l, e_r, lr;
        bool undecided;

        Point m = midpoint(x->max[lower].max(),y->min[lower].min());

        while (!(isLeaf(x) && isLeaf(y))) {
            undecided = true;
            e_l = x->val[lower];
            e_r = y->val[lower];
            lr = Bridge(midpoint(e_l),midpoint(e_r));
            if (!isLeaf(x) && slope_comp<lower>(e_l,lr)){
                x = stepLeft<lower>(x); undecided = false;
            }
            if (!isLeaf(y) && slope_comp<lower>(lr,e_r)) {
                y = stepRight<lower>(y); undecided = false;
            }
            if (undecided) {
                if (!isLeaf(x) && m_comp<lower>(e_l,e_r,m) || isLeaf(y)) {
                    x = stepRight<lower>(x);
                } else {
                    y = stepLeft<lower>(y);
                }
            }
        }
        return Bridge(x->val[lower].min(),y->val[lower].max());
    }

    void onUpdate(Node* x){
        if(isLeaf(x)) return;
        x->val[0] = findBridge<false>(x);
        x->val[1] = findBridge<true>(x);
    }

    template<bool lower>
    inline
    Node* stepLeft(Node* v){
        auto x = v->val[lower].min().x();
        v = v->left;
        while(v && v->val[lower].max().x() > x) v = v->left;
        return v;
    }

    template<bool lower>
    inline
    Node* stepRight(Node* v){
        auto x = v->val[lower].max().x();
        v = v->right;
        while(v && v->val[lower].min().x() < x) v = v->right;
        return v;
    }

    template<bool lower>
    Node* find(const Point key, const bool left){
        auto current = AVLTree<Bridges<Traits>>::root;
        while(current && !isLeaf(current)){
            if(current->val[lower][!left] == key) return current;
            else if (current->val[lower].min().x() < key.x()) current = current->right;
            else current = current->left;
        }
        return nullptr;
    }

    template<bool lower>
    Node* hullSuccessor(const Point key){
        auto current = AVLTree<Bridges<Traits>>::root;
        Point min;
        while(current){
            min = current->val[lower].min();
            if(min == key) break;
            else if (min < key) current = stepRight<lower>(current);
            else current = stepLeft<lower>(current);
        }
        return current;
    }

    template<bool lower>
    Node* hullSuccessor(const Node* p){
        return hullSuccessor<lower>(p->val[lower].max());
    }

    template<bool lower>
    Node* hullPredecessor(const Point key){
        auto current = AVLTree<Bridges<Traits>>::root;
        Point max;
        while(current){
            max = current->val[lower].max();
            if(max == key) break;
            else if (max < key) current = stepRight<lower>(current);
            else current = stepLeft<lower>(current);
        }
        return current;
    }

    template<bool lower>
    Node* hullPredecessor(const Node* p){
        return hullPredecessor<lower>(p->val[lower].min());
    }

    template<bool lower>
    bool covers(Point p){
        auto current = AVLTree<Bridges<Traits>>::root;
        while(current){
            if(current->val[lower].min() <= p){
                if(p <= current->val[lower].max()){
                    return cover_comp<lower>(p,current->val[lower]);
                } else current = stepRight<lower>(current);
            } else current = stepLeft<lower>(current);
        }
        return false;
    }

    template<bool lower>
    std::vector<Point> hullPoints(){
        std::vector<Point> res;
        Node* e = AVLTree<Bridges<Traits>>::root;
        if(!e || isLeaf(e)) return res;
        while(e && !isLeaf(e)){
            res.insert(res.begin(), e->val[lower].min());
            e = hullPredecessor<lower>(res.front());
        }
        e = AVLTree<Bridges<Traits>>::root;
        while(e && !isLeaf(e)){
            res.push_back(e->val[lower].max());
            e = hullSuccessor<lower>(res.back());
        }
        if(lower) res.pop_back();
        else {
            res.erase(res.begin());
            std::reverse(res.begin(), res.end());
        }
        return res;
    }

public:
    void insert(Point p){
        AVLTree<Bridges<Traits>>::insert({Bridge(p,p),Bridge(p,p)});
    }

    void remove(Point p){
        AVLTree<Bridges<Traits>>::remove({Bridge(p,p),Bridge(p,p)});
    }

    bool covers(Point p){
        return covers<false>(p) && covers<true>(p);
    }

    std::vector<Point> upperHullPoints(){
        return hullPoints<false>();
    }

    std::vector<Point> lowerHullPoints(){
        return hullPoints<true>();
    }
};
#endif //DYNAMICCONVEXHULL_CHTREE_H
