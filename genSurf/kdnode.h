#ifndef KDNODE_H
#define KDNODE_H

#include "reconstruct.h"

template <class T> class KDNode {
  public:
    KDNode();
    KDNode(const T& item, KDNode<T>* nextPtr = NULL);
    T data;
    KDNode<T>* NextKDNode(); // for next node
    void InsertAfter(KDNode<T>* p);
    KDNode<T>* DeleteAfter();
    KDNode<T>* GetNode(const T& item, KDNode<T>* nextPtr = NULL);

  private:
    KDNode<T>* next;
};

#endif // KDNODE_H