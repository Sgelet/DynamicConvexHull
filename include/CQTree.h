//
// Created by etoga on 5/12/23.
//

#ifndef DYNAMICCONVEXHULL_CQTREE_H
#define DYNAMICCONVEXHULL_CQTREE_H

#include <vector>
#include "AvlTree.h"
#include "util.h"
#define isLeaf this->isLeaf

template<class Traits>
struct comp_xy{
    static constexpr auto comp = typename Traits::Compare_xy_2();
    bool operator()(const typename Traits::Segment_2& lhs, const typename Traits::Segment_2& rhs) const {
        bool a = comp(lhs[0],rhs[0]) == CGAL::SMALLER;
        bool b = comp(lhs[1],rhs[1]) == CGAL::EQUAL;
        return comp(lhs[0],rhs[0]) == CGAL::SMALLER || (comp(lhs[1],rhs[1]) == CGAL::EQUAL && comp(lhs[0],rhs[0]) != CGAL::EQUAL);
    }
};

template<class Traits>
using CQueue = AVLTree<typename Traits::Segment_2,comp_xy<Traits>>;

template<class Traits>
struct BCQ{
    using Bridge = typename Traits::Segment_2;
    Bridges<Traits> bridges;
    std::array<CQueue<Traits>,2> hulls;

    BCQ(Bridge x, Bridge y):bridges(x,y){};

    BCQ(BCQ& b):bridges(b.bridges){};

    constexpr bool operator==(const BCQ& b) const{
        return bridges == b.bridges;
    }

    constexpr bool operator!=(const BCQ& b) const{
        return !(*this == b);
    }

    constexpr bool operator<(const BCQ& b) const{
        return bridges < b.bridges;
    }

    constexpr bool operator<=(const BCQ& b) const{
        return !(b < *this);
    }
};

template<class Traits>
class CQTree : AVLTree<BCQ<Traits>>{
using Node = typename AVLTree<BCQ<Traits>>::Node;
using HNode = typename CQueue<Traits>::Node;
using Bridge = typename Traits::Segment_2;
using Point = typename Traits::Point_2;
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
    bool m_comp(const Bridge& l, const Bridge& r, const Point& m) {
        if(l.is_vertical()) return lower;
        if(r.is_vertical()) return !lower;
        CGAL::Comparison_result res = compare_at_x(m, l.supporting_line(), r.supporting_line());
        if (res == CGAL::EQUAL) return true;
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
        HNode* x = v->left->val.hulls[lower].root;
        HNode* y = v->right->val.hulls[lower].root;
        Bridge e_l, e_r, lr;
        bool undecided;
        bool foundl = false;
        bool foundr = false;
        Point l,r;
        Point m = midpoint(v->left->max.bridges[lower].max(),v->right->min.bridges[lower].min());
        if(!x){
            foundl = true;
            l = v->left->val.bridges[lower].min();
        }
        if(!y){
            foundr = true;
            r = v->right->val.bridges[lower].min();
        }
        while (!foundl || !foundr) {
            undecided = true;
            if(!foundl){
                e_l = x->val;
                l = midpoint(x->val);
            }
            if(!foundr) {
                e_r = y->val;
                r = midpoint(y->val);
            }
            lr = Bridge(l,r);
            if (!foundl && slope_comp<lower>(e_l,lr)){
                if(x->left) x = x->left;
                else{
                    foundl = true;
                    l = x->val.min();
                }
                undecided = false;
            }
            if (!foundr && slope_comp<lower>(lr,e_r)) {
                if(y->right) y = y->right;
                else{
                    foundr = true;
                    r = y->val.max();
                }
                undecided = false;
            }
            if (undecided) {
                if (foundr ||  (!foundl && m_comp<lower>(e_l,e_r,m))) {
                    if(x->right) x = x->right;
                    else{
                        foundl = true;
                        l = x->val.max();
                    }
                } else {
                    if(y->left) y = y->left;
                    else{
                        foundr = true;
                        r = y->val.min();
                    }
                }
            }
        }
        return {l,r};
    }

    void onUpdate(Node* x){
        if(isLeaf(x)) return;
        x->val.bridges[0] = findBridge<false>(x);
        x->val.bridges[1] = findBridge<true>(x);
        CQueue<Traits> left, right;
        for(int i=0; i<2; ++i){
            x->left->val.hulls[i].split({x->val.bridges[i][0],x->val.bridges[i][0]},&left);
            x->right->val.hulls[i].split({x->val.bridges[i][1],x->val.bridges[i][1]},&right);
            x->left->val.hulls[i].join(x->val.bridges[i],&right);
            x->val.hulls[i].join(&(x->left->val.hulls[i]));
            x->left->val.hulls[i].join(&left);
        }
    }

    void onVisit(Node* x){
       if(isLeaf(x)) return;
       CQueue<Traits> right;
       for(int i=0;i<2;++i){
           x->val.hulls[i].split(x->val.bridges[i],&right);
           x->val.hulls[i].join(&(x->left->val.hulls[i]));
           x->left->val.hulls[i].join(&(x->val.hulls[i]));
           x->right->val.hulls[i].join(&right);
       }
    }

    void hullPoints(HNode* x, std::vector<Point>& acc){
        if(!x) return;
        hullPoints(x->left,acc);
        acc.push_back(x->val[0]);
        hullPoints(x->right,acc);
    }

    void hullPoints2(HNode* x, std::vector<Point>& acc){
       if(!x) return;
       hullPoints2(x->right,acc);
       acc.push_back(x->val[1]);
       hullPoints2(x->left,acc);
    }

    template<bool lower>
    bool covers(Point p){
        auto current = AVLTree<BCQ<Traits>>::root->val.hulls[lower].root;
        while(current){
            if(current->val.min().x() <= p.x()){
                if(p.x() <= current->val.max().x()){
                    return cover_comp<lower>(p,current->val);
                } else current = current->right;
            } else current = current->left;
        }
        return false;
    }

public:
    void insert(Point p){
        AVLTree<BCQ<Traits>>::insert(BCQ<Traits>(Bridge(p,p),Bridge(p,p)));
    }

    void remove(Point p){
        AVLTree<BCQ<Traits>>::remove(BCQ<Traits>(Bridge(p,p),Bridge(p,p)));
    }

    bool covers(Point p){
        return covers<false>(p) && covers<true>(p);
    }

    std::vector<Point> upperHullPoints(){
        std::vector<Point> res;
        if(this->root) hullPoints2(this->root->val.hulls[0].root,res);
        return res;
    }
    std::vector<Point> lowerHullPoints(){
        std::vector<Point> res;
        if(this->root) hullPoints(this->root->val.hulls[1].root,res);
        return res;
    }
};
#endif //DYNAMICCONVEXHULL_CQTREE_H
