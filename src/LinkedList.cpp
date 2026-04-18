#include "LinkedList.h"

LinkedList::LinkedList()
    : head_(nullptr) {}

// Constructeur
LinkedList::LinkedList(LinkedList&& other) noexcept
    : head_(other.head_) {
    other.head_ = nullptr;
}

LinkedList& LinkedList::operator=(LinkedList&& other) noexcept {
    if (this != &other) {

        Node* cur = head_;
        while (cur != nullptr) {
            Node* nxt = cur->next;
            delete cur;
            cur = nxt;
        }

        head_ = other.head_;
        other.head_ = nullptr;
    }
    return *this;
}

// Destructeur
LinkedList::~LinkedList() {
    Node* cur = head_;
    while (cur != nullptr) {
        Node* nxt = cur->next;
        delete cur;
        cur = nxt;
    }
}

void LinkedList::insert(uint64_t kmerCode, const GenomicPosition& position) {
    Node* n = new Node(kmerCode, position);
    n->next = head_;
    head_ = n;
}

Node* LinkedList::first() const {
    return head_;
}
